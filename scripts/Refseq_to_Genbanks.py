import pandas as pd
from Bio import Entrez
import subprocess
import sys

def search_nucleotide(refseq_accession):
    Entrez.email = ""  
    handle = Entrez.esearch(db="nucleotide", term=refseq_accession, idtype="acc")
    record = Entrez.read(handle)
    handle.close()
    if record['Count'] and int(record['Count']) > 0:
        nucleotide_id = record['IdList'][0]
        return nucleotide_id
    else:
        return None

def download_genbank(nucleotide_id):
    Entrez.email = ""  
    handle = Entrez.efetch(db="nucleotide", id=nucleotide_id, rettype="gb", retmode="text")
    genbank_data = handle.read()
    handle.close()
    return genbank_data

def main(excel_file):
    # Read Excel input file
    df = pd.read_excel(excel_file)

    # Iterate over each RefSeq accession, search and download GenBank file, and concatenate them
    extracted_genbanks = ""
    for refseq_accession in df['RefSeq Accessions']:
        nucleotide_id = search_nucleotide(refseq_accession)
        if nucleotide_id:
            genbank_data = download_genbank(nucleotide_id)
            extracted_genbanks += genbank_data

    # Write concatenated GenBank data to output file
    output_file = "extracted_genbanks.gb"
    with open(output_file, "w") as f:
        f.write(extracted_genbanks)

    print("GenBank files extracted and concatenated to", output_file)
    
    # Call the next script and pass the necessary files as arguments
    next_script = "3UTR_extraction.py"
    subprocess.run(["python", next_script, output_file, excel_file])

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python first_script.py <path_to_excel_file>")
        sys.exit(1)

    excel_file = sys.argv[1]
    main(excel_file)
