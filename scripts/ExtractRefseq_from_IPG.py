import pandas as pd
from Bio import Entrez
import subprocess 

def search_ipg(accession_number):
    Entrez.email = ""  

    # Search NCBI Identical Protein Groups (IPG) database for the accession number
    handle = Entrez.esearch(db="protein", term=f"{accession_number}[Accession]", retmax=1)
    record = Entrez.read(handle)
    handle.close()

    if "Count" in record and int(record["Count"]) > 0:
        ipg_id = record["IdList"][0]
        # Fetch IPG summary for the accession number
        handle = Entrez.efetch(db="protein", id=ipg_id, rettype="ipg", retmode="xml")
        ipg_record = Entrez.read(handle)
        handle.close()

        return ipg_record
    else:
        return None

def extract_ipg_info(accession_number, ipg_record):
    ipg_info = {"Protein Accession": accession_number}
    refseq_accession = None

    if ipg_record and 'IPGReport' in ipg_record:
        # Search for 'ProteinList' and extract the first 'accver' following it
        if 'ProteinList' in ipg_record['IPGReport']:
            for protein in ipg_record['IPGReport']['ProteinList']:
                if 'CDSList' in protein:
                    for cds in protein['CDSList']:
                        if 'accver' in cds.attributes:
                            refseq_accession = cds.attributes['accver']
                            break  # We only need the first one
                if refseq_accession:
                    break  # Stop once we've found the RefSeq Accession

    ipg_info["RefSeq Accessions"] = refseq_accession

    return ipg_info

def main(input_file):
    # Read input Excel file
    df = pd.read_excel(input_file)

    ipg_info_list = []
    # Iterate over protein accession numbers and search NCBI IPG
    for accession_number in df["protein accession"]:
        print(f"Searching IPG for accession number: {accession_number}")
        ipg_record = search_ipg(accession_number)
        if ipg_record:
            ipg_info = extract_ipg_info(accession_number, ipg_record)
            if ipg_info:
                ipg_info_list.append(ipg_info)

    # Create DataFrame from IPG info list
    ipg_df = pd.DataFrame(ipg_info_list)

    # Write IPG info to output Excel file
    output_file = "refseqs.xlsx"
    ipg_df.to_excel(output_file, index=False)

    # Call the next script and pass the output Excel file as an argument
    next_script = "Refseq_to_Genbanks.py"
    subprocess.run(["python", next_script, output_file])

if __name__ == "__main__":
    input_file = input("Please provide the input Excel file name (or path): ")
    main(input_file)
