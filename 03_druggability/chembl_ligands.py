#!/usr/bin/env python3
"""Retrieve annotated small-molecule ligands for the hub genes from ChEMBL.

For each hub gene the script resolves human SINGLE PROTEIN targets, keeps the
IC50 bioactivities that carry a pChEMBL value, and deduplicates to one row per
molecule per target. The result is the evidence behind the target selection:
which of the prioritised genes actually have a tractable chemistry, and how the
compound taken into the experiment compares with the rest.

Two details matter and were wrong in the earlier version of this step.

*Targets are matched on an exact gene symbol, not a substring.* ChEMBL's
``target_synonym__icontains`` is a substring search, so a query for HDAC1 also
returns HDAC10 and HDAC11 and inflates the count for HDAC1 roughly threefold.
Each candidate target is therefore kept only if one of its component synonyms
equals the gene symbol exactly.

*Bioactivities are restricted to IC50 with a pChEMBL value.* Mixing Ki, EC50
and IC50 makes potencies incomparable across targets, and records without a
pChEMBL value have no standardised potency at all.

ChEMBL is a live database, so the counts reflect the day the query ran; the
retrieval date is written into the output.

Usage:  python3 03_druggability/chembl_ligands.py
"""

import datetime
import os
import statistics

import pandas as pd
from chembl_webresource_client.new_client import new_client
from chembl_webresource_client.settings import Settings

Settings.Instance().TIMEOUT = 20

HUBS = "results/ppi/intersected_genes_up_alz_hub.txt"
HUBS_DOWN = "results/ppi/intersected_genes_down_alz_hub.txt"
OUT = "results/druggability"

def hub_genes():
    genes = []
    for path in (HUBS, HUBS_DOWN):
        with open(path) as fh:
            genes += [l.strip() for l in fh if l.strip() and l.strip() != "gene"]
    return genes


def exact_targets(gene):
    """Human SINGLE PROTEIN targets whose official gene symbol is `gene`.

    Two filters are needed and dropping either one lets the wrong protein in.
    The symbol must match exactly, because ``target_synonym__icontains`` is a
    substring search that returns HDAC10 and HDAC11 for a query of HDAC1. The
    matching synonym must also be typed ``GENE_SYMBOL`` rather than
    ``GENE_SYMBOL_OTHER``: ChEMBL lists PRKACA among the aliases of protein
    kinase C alpha (CHEMBL299, UniProt P17252), a different protein from the
    cAMP-dependent kinase the symbol actually names (CHEMBL4101, P17612), and
    an untyped match assigns that target's ligands to PRKACA.
    """
    hits = new_client.target.filter(target_synonym__icontains=gene,
                                    organism="Homo sapiens",
                                    target_type="SINGLE PROTEIN")
    out = []
    for t in hits:
        official = {s.get("component_synonym", "").upper()
                    for comp in t.get("target_components", [])
                    for s in comp.get("target_component_synonyms", [])
                    if s.get("syn_type") == "GENE_SYMBOL"}
        if gene.upper() in official:
            acc = next((c.get("accession", "") for c in t.get("target_components", [])), "")
            out.append((t["target_chembl_id"], t["pref_name"], acc))
    return out


def main():
    os.makedirs(OUT, exist_ok=True)
    today = datetime.date.today().isoformat()
    rows = []
    for gene in hub_genes():
        targets = exact_targets(gene)
        if not targets:
            rows.append(dict(gene=gene, target_chembl_id="", target_pref_name="",
                             uniprot="", molecule_chembl_id="",
                             pchembl_value="", n_measurements=""))
            print("  {:<9} no human single-protein target".format(gene))
            continue
        for tid, pref, acc in targets:
            acts = new_client.activity.filter(
                target_chembl_id=tid, standard_type="IC50",
                pchembl_value__isnull=False).only(
                ["molecule_chembl_id", "pchembl_value"])
            # A molecule is usually measured many times against the same
            # target and the values spread widely (vorinostat against HDAC1
            # spans pChEMBL 4.5 to 9.6 across 416 records), so each molecule is
            # summarised by the median of its measurements rather than by
            # whichever record the API happened to return first.
            seen = {}
            for a in acts:
                mid, val = a["molecule_chembl_id"], a["pchembl_value"]
                if mid and val is not None:
                    seen.setdefault(mid, []).append(float(val))
            if not seen:
                # A target with no qualifying bioactivity still belongs in the
                # table: "queried, nothing found" and "never queried" are
                # different claims, and without this row the gene disappears.
                rows.append(dict(gene=gene, target_chembl_id=tid,
                                 target_pref_name=pref, uniprot=acc,
                                 molecule_chembl_id="", pchembl_value="",
                                 n_measurements=""))
            for mid, vals in seen.items():
                rows.append(dict(gene=gene, target_chembl_id=tid,
                                 target_pref_name=pref, uniprot=acc,
                                 molecule_chembl_id=mid,
                                 pchembl_value=statistics.median(vals),
                                 n_measurements=len(vals)))
            print("  {:<9} {:<14} {:<38} {} ligands".format(
                gene, tid, pref[:36], len(seen)))

    df = pd.DataFrame(rows)
    df["retrieved"] = today
    df.to_csv(os.path.join(OUT, "hub_gene_ligands.tsv"), sep="\t", index=False)

    counts = (df[df.molecule_chembl_id != ""]
              .groupby("gene").molecule_chembl_id.nunique()
              .reindex(hub_genes(), fill_value=0).sort_values(ascending=False))
    counts.to_csv(os.path.join(OUT, "hub_gene_ligand_counts.tsv"), sep="\t",
                  header=["n_ligands"])
    print("\nwrote {}/hub_gene_ligands.tsv  ({} rows, retrieved {})".format(
        OUT, len(df), today))
    print("wrote {}/hub_gene_ligand_counts.tsv".format(OUT))


if __name__ == "__main__":
    main()
