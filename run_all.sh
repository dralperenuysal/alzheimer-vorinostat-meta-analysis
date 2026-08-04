#!/usr/bin/env bash
#
# Run the published pipeline end to end, or one stage of it.
#
#     bash run_all.sh              # list the stages and stop
#     bash run_all.sh all          # everything, start to finish
#     bash run_all.sh 01 02       # just the meta-analysis and the PPI step
#     bash run_all.sh from 02      # stage 02 onwards
#
# The chain runs from prepared expression matrices to the hub genes and on to
# the qPCR analysis and figures. Primer design sits outside it: see the comment
# above the stage list and 05_primer_design/README.md.
#
# Stage 01 reads data/, which is not committed. Prepare it first; the README
# gives the expected file names and layout.
#
# Everything runs in minutes. Stages 02 and 03 need network access to STRING
# and ChEMBL.
#
# After a run, compare what came out against the published results:
#     python3 verify_reproduction.py

set -uo pipefail

cd "$(dirname "$0")"

RED=$'\033[31m'; GREEN=$'\033[32m'; BOLD=$'\033[1m'; DIM=$'\033[2m'; OFF=$'\033[0m'
[ -t 1 ] || { RED=""; GREEN=""; BOLD=""; DIM=""; OFF=""; }

STAGES=(01 02 03 04)

declare -A STAGE_NAME=(
    [01]="meta-analysis: DExMA random effects, then the intersection"
    [02]="PPI network and hub genes (STRING, needs network access)"
    [03]="druggability: ChEMBL ligands for the hub genes (needs network access)"
    [04]="figures: volcano plots, Venn diagrams, PPI networks, druggability, qPCR panels"
)

# 05_primer_design/ is deliberately not a stage. Primer design is bench work
# that happened once, between the hub genes and the qPCR, and its result is the
# primer table in the manuscript. It is not part of the computational chain and
# is not reproducible in the sense the rest of this script means: the Primer3
# build behind the published primers was not recorded, and four of the fourteen
# genes were designed in NCBI's web interface rather than by the script. The
# templates are pinned by accession, so what each published pair was designed on
# is known exactly; the design step itself is not re-run.
# See 05_primer_design/README.md.

FAILED=()
STARTED=$(date +%s)

announce() { printf '\n%s== stage %s: %s%s\n' "$BOLD" "$1" "${STAGE_NAME[$1]}" "$OFF"; }

# Run one command, report it, and keep going if it fails. The summary at the
# end lists everything that went wrong.
step() {
    printf '%s-- %s%s\n' "$DIM" "$*" "$OFF"
    if ! "$@"; then
        printf '%s   FAILED: %s%s\n' "$RED" "$*" "$OFF"
        FAILED+=("$*")
    fi
}

stage_01() {
    announce 01
    step Rscript 01_meta_analysis/DExMA_meta_microarray.R
    step Rscript 01_meta_analysis/DExMA_meta_rnaseq.R
    step Rscript 01_meta_analysis/result_intersection.R
}

stage_02() {
    announce 02
    local direction
    for direction in up down; do
        step python3 02_ppi_network/string_network.py \
            --genes "results/meta_analysis/intersected_genes_${direction}_alz.txt" \
            --out results/ppi
        step python3 02_ppi_network/hub_genes.py \
            --edges "results/ppi/intersected_genes_${direction}_alz_edges.tsv" \
            --top 14 --out results/ppi
    done
}

stage_03() {
    announce 03
    step python3 03_druggability/chembl_ligands.py
}

stage_04() {
    announce 04
    step Rscript 04_figures/volcano_plots.R
    step Rscript 04_figures/venn_diagrams.R
    step python3 04_figures/make_network_figures.py
    step python3 04_figures/make_druggability_figure.py
    step python3 04_figures/make_qpcr_figures.py
}

usage() {
    echo "Usage: bash run_all.sh [all | from <stage> | <stage> ...]"
    echo
    echo "Stages:"
    local stage
    for stage in "${STAGES[@]}"; do
        printf '  %-4s %s\n' "$stage" "${STAGE_NAME[$stage]}"
    done
    echo
    echo "Stage 01 reads prepared matrices from data/; see the README for the"
    echo "expected layout. Stages 02 and 03 need network access."
    echo
    echo "Primer design (05_primer_design/) is not a stage. It was a one-off"
    echo "bench step between the hub genes and the qPCR; see its README."
}

# ------------------------------------------------------------ stage selection

[ $# -eq 0 ] && { usage; exit 0; }

SELECTED=()
case "$1" in
    all)
        SELECTED=("${STAGES[@]}")
        ;;
    from)
        [ $# -ge 2 ] || { echo "from what? e.g. 'bash run_all.sh from 01'" >&2; exit 1; }
        collecting=0
        for stage in "${STAGES[@]}"; do
            [ "$stage" = "$2" ] && collecting=1
            [ $collecting -eq 1 ] && SELECTED+=("$stage")
        done
        [ ${#SELECTED[@]} -gt 0 ] || { echo "Unknown stage '$2'." >&2; usage >&2; exit 1; }
        ;;
    -h|--help|help)
        usage; exit 0
        ;;
    *)
        for requested in "$@"; do
            found=0
            for stage in "${STAGES[@]}"; do
                [ "$stage" = "$requested" ] && { SELECTED+=("$stage"); found=1; }
            done
            [ $found -eq 1 ] || { echo "Unknown stage '$requested'." >&2; usage >&2; exit 1; }
        done
        ;;
esac

for command in Rscript python3; do
    command -v "$command" >/dev/null || {
        echo "${RED}$command is not on PATH. Create the environments in env/ first.${OFF}" >&2
        exit 1
    }
done

echo "${BOLD}Running stages: ${SELECTED[*]}${OFF}"

for stage in "${SELECTED[@]}"; do
    "stage_$stage"
done

# ------------------------------------------------------------------- summary

ELAPSED=$(( $(date +%s) - STARTED ))
printf '\n%s%s%s\n' "$BOLD" "$(printf '=%.0s' {1..72})" "$OFF"
printf 'Finished stages %s in %dm %ds\n' "${SELECTED[*]}" $((ELAPSED / 60)) $((ELAPSED % 60))

if [ ${#FAILED[@]} -eq 0 ]; then
    printf '%sEvery step completed.%s\n' "$GREEN" "$OFF"
else
    printf '%s%d step(s) failed:%s\n' "$RED" "${#FAILED[@]}" "$OFF"
    for failure in "${FAILED[@]}"; do
        printf '  %s\n' "$failure"
    done
    printf '\nLater stages read what earlier ones wrote, so a failure here usually\n'
    printf 'means the results below it are incomplete rather than wrong.\n'
fi

printf '\nNow compare against the published results:\n  python3 verify_reproduction.py\n'
[ ${#FAILED[@]} -eq 0 ]
