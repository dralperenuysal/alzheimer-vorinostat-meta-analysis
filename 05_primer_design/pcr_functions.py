def design_qpcr_primers(sequence, target_start=0, target_length=150):
    """
    Designs qPCR-compatible primers for a target sequence in the HDAC1 gene.

    Parameters:
    - sequence (str): DNA sequence of the target region in the HDAC1 gene.
    - target_start (int): Start position of the target region within the sequence.
    - target_length (int): Length of the target region to amplify.

    Returns:
    - dict: Primer design output with primer sequences and additional metadata.
    """

    from primer3 import design_primers

    # Primer3 settings for qPCR-compatible primer design
    settings = {
        'PRIMER_PRODUCT_SIZE_RANGE': [[50, 150]],  # Optimal range for qPCR
        'PRIMER_OPT_SIZE': 20,
        'PRIMER_MIN_SIZE': 18,
        'PRIMER_MAX_SIZE': 24,
        'PRIMER_MIN_TM': 57.0,
        'PRIMER_OPT_TM': 60.0,
        'PRIMER_MAX_TM': 63.0,
        'PRIMER_MIN_GC': 40.0,
        'PRIMER_OPT_GC_PERCENT': 50.0,
        'PRIMER_MAX_GC': 60.0,
        'PRIMER_MAX_POLY_X': 5,
        'PRIMER_MAX_SELF_ANY': 8.0,
        'PRIMER_MAX_SELF_END': 3.0,
        'PRIMER_SALT_CONC': 50.0,
        'PRIMER_DNA_CONC': 200.0,
    }

    # Primer design with Primer3
    result = design_primers(
        {
            'SEQUENCE_ID': 'HDAC1_target',
            'SEQUENCE_TEMPLATE': sequence,
            'SEQUENCE_TARGET': [target_start, target_length],
        },
        settings
    )

    # Output the result
    if result:
        print(f"Left Primer: {result['PRIMER_LEFT_0_SEQUENCE']}")
        print(f"Right Primer: {result['PRIMER_RIGHT_0_SEQUENCE']}")
        print(f"Product Size: {result['PRIMER_PAIR_0_PRODUCT_SIZE']}")
        print(f"Left Primer Tm: {result['PRIMER_LEFT_0_TM']}")
        print(f"Right Primer Tm: {result['PRIMER_RIGHT_0_TM']}")
        print(f"Left Primer GC Content: {result['PRIMER_LEFT_0_GC_PERCENT']}")
        print(f"Right Primer GC Content: {result['PRIMER_RIGHT_0_GC_PERCENT']}")
        return result
    else:
        print("No suitable primers found with given parameters.")
        return None

def download_fasta(gene, organism, email, saving_loc):
    import os
    from Bio import Entrez, SeqIO
    # Set the email address for NCBI (required by NCBI usage policy)
    Entrez.email = email  # Replace with your email

    # Define gene name and organism
    gene_name = gene  # Example gene: HDAC1
    organism = organism  # Example organism: Homo sapiens (human)
    file = os.path.join(saving_loc, f"{gene_name}_genomic.fasta") # define the saving file

    # Step 1: Search for the gene's ID in the NCBI nucleotide database
    search_handle = Entrez.esearch(db="nucleotide", term=f"{gene_name}[Gene] AND {organism}[Organism]", rettype="fasta")
    search_results = Entrez.read(search_handle)
    search_handle.close()

    # Step 2: Fetch the genomic sequence using the first result ID
    if search_results["IdList"]:
        gene_id = search_results["IdList"][0]  # Select the first result ID
        fetch_handle = Entrez.efetch(db="nucleotide", id=gene_id, rettype="fasta", retmode="text")
        gene_sequence = fetch_handle.read()
        fetch_handle.close()

        # Step 3: Save the sequence as a FASTA file
        with open(file, "w") as fasta_file:
            fasta_file.write(gene_sequence)
        print(f"{gene_name}_genomic.fasta file has been saved.")
    else:
        print("Gene not found.")