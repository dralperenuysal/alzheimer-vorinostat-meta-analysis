from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import os

files = os.listdir("data/qPCR")

for file in files:
    # Primer dosyasını aç ve sırayı oku
    f = os.path.join("data/qPCR", file)
    with open(f) as fasta_file:
        fasta_data = fasta_file.read()

    # NCBI BLAST web API'yi kullanarak BLAST çalıştır
    result_handle = NCBIWWW.qblast("blastn", "nt", fasta_data)

    # Sonuçları kaydet
    name = file.split("_")[0]
    save = f"results/qpcr/{name}_blast_results.xml"
    with open(save, "w") as out_file:
        out_file.write(result_handle.read())

    print(f"BLAST results saved to {save}")

results = os.listdir("data/qPCR")
results = [res for res in results if res.endswith(".xml")]

for res in results:
    from Bio.Blast import NCBIXML

    # XML dosyasını aç ve parse et
    with open(f"results/qpcr/{res}") as result_handle:
        blast_records = NCBIXML.parse(result_handle)

        # XML içindeki her bir BLAST kaydını dolaş
        for blast_record in blast_records:
            print(f"Query ID: {blast_record.query_id}")

            for alignment in blast_record.alignments:
                print(f"  Subject ID: {alignment.hit_id}")
                for hsp in alignment.hsps:
                    print(f"    E-value: {hsp.expect}")
                    print(f"    Score: {hsp.score}")
                    print(f"    Identities: {hsp.identities}")
                    print(f"    Query Start: {hsp.query_start}")
                    print(f"    Query End: {hsp.query_end}")
                    print(f"    Subject Start: {hsp.sbjct_start}")
                    print(f"    Subject End: {hsp.sbjct_end}")
                    print()

import primer3

primer3.design_primers()
