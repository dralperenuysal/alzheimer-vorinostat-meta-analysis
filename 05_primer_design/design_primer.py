"""Earlier, incomplete pass at primer design. NOT the route used in the study.

The primers reported in the manuscript were designed with design_primer.R.
This script differs from it in ways that make the two non-interchangeable:

  - it asks Primer3 for 70-150 bp amplicons, where the study used 150-200
  - it skips PSMD14 as too long
  - it never writes its results anywhere, so it cannot have produced
    results/qpcr/*_primers.tsv

It is kept because it is part of how the analysis was arrived at, not because
it is a working alternative. Use design_primer.R.
"""

import os
import primer3
from Bio import Entrez, SeqIO
import pandas as pd
from pcr_functions import download_fasta
from primer3 import design_primers

# read the upregulated genes
genes_all = dict()
up_loc = "results/ppi/intersected_genes_up_alz_hub.txt"
down_loc = "results/ppi/intersected_genes_down_alz_hub.txt"
genes_all["up_genes"] = pd.read_csv(up_loc, header="infer").gene.tolist()
genes_all["down_genes"] = pd.read_csv(down_loc, header="infer").gene.tolist()

# download the fasta files with a for loop
# NCBI Entrez requires a contact address with every request.
ENTREZ_EMAIL = os.environ.get("ENTREZ_EMAIL", "your.name@example.com")

for pattern in genes_all.keys():
    print(pattern)
    for gene in genes_all[pattern]:
        print(f"{gene} is being downloaded.")
        download_fasta(gene=gene,
                       organism="Homo sapiens",
                       email=ENTREZ_EMAIL,
                       saving_loc="data/raw/qPCR")
        print("-----------")

# design the primers with a for loop
for pattern in genes_all.keys():
    print(pattern)
    for gene in genes_all[pattern]:
        if gene != "PSMD14": # PSMD14 is too long
            # read the sequence from fasta file
            fasta_file = f"data/raw/qPCR/{gene}_genomic.fasta"
            record = SeqIO.read(fasta_file, "fasta")
            sequence = str(record.seq)

            # Sequence-specific arguments
            seq_args = {
                'SEQUENCE_ID': 'HDAC1_gene',
                'SEQUENCE_TEMPLATE': sequence,
                'SEQUENCE_INCLUDED_REGION': [0, len(sequence)],  # Target the full sequence
            }

            # Global arguments for qPCR primer design
            global_args = {
                'PRIMER_PRODUCT_SIZE_RANGE': [[70, 150]],      # Suitable range for qPCR products
                'PRIMER_OPT_SIZE': 20,
                'PRIMER_MIN_SIZE': 18,
                'PRIMER_MAX_SIZE': 24,
                'PRIMER_OPT_TM': 60.0,
                'PRIMER_MIN_TM': 57.0,
                'PRIMER_MAX_TM': 63.0,
                'PRIMER_MIN_GC': 40.0,
                'PRIMER_OPT_GC_PERCENT': 50.0,
                'PRIMER_MAX_GC': 60.0,
                'PRIMER_MAX_POLY_X': 5,                        # Max mononucleotide repeat length
                'PRIMER_MAX_SELF_ANY': 8.0,                    # Max allowable primer self-complementarity
                'PRIMER_MAX_SELF_END': 3.0,
                'PRIMER_SALT_CONC': 50.0,
                'PRIMER_DNA_CONC': 200.0,
            }

            # Optional mispriming and mishybridization libraries
            misprime_lib = {
                'human_misprime': 'AGCTTGACCTGACCCAGTGAAGCTGTGAACTTCCAGAACGCGGAGGAGG'
            }
            mishyb_lib = {
                'human_mishyb': 'GATGCTGACCTGAATGGCCTGAGCTGAAGGAATGGAAGAGCTGAAG'
            }

            # Run primer design
            result = design_primers(seq_args, global_args, misprime_lib, mishyb_lib)

            # Output the primer design results
            print(f"Results for {gene}")
            print("Left Primer:", result['PRIMER_LEFT_0_SEQUENCE'])
            print("Right Primer:", result['PRIMER_RIGHT_0_SEQUENCE'])
            print("Product Size:", result['PRIMER_PAIR_0_PRODUCT_SIZE'])
            print("Left Primer Tm:", result['PRIMER_LEFT_0_TM'])
            print("Right Primer Tm:", result['PRIMER_RIGHT_0_TM'])
            print("GC content of Left Primer:", result['PRIMER_LEFT_0_GC_PERCENT'])
            print("GC content of Right Primer", result['PRIMER_RIGHT_0_GC_PERCENT'])
            print("\n", "-----")

