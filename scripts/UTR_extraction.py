from Bio import SeqIO
import pandas as pd
import sys 

def extract_upstream_downstream_sequence(genbank_file, protein_ids_file):
    # Read the protein IDs from the Excel file
    protein_ids_df = pd.read_excel(protein_ids_file)

    # Check if the required columns exist in the Excel file
    if 'Protein Accession' not in protein_ids_df.columns or 'RefSeq Accessions' not in protein_ids_df.columns:
        print("The required columns are missing from the Excel file.")
        print("Available columns are:", protein_ids_df.columns)
        return
    
    # Create a mapping from Protein Accession to RefSeq Accessions
    accession_map = dict(zip(protein_ids_df['Protein Accession'], protein_ids_df['RefSeq Accessions']))
    
    print("Protein IDs and their corresponding RefSeq Accessions:", accession_map)

    # Prepare lists to store data for Excel output
    data = []

    # Open the output FASTA files
    with open("5UTR.fasta", "w") as upstream_fasta, open("3UTR.fasta", "w") as downstream_fasta:
        # Iterate through the GenBank file
        for record in SeqIO.parse(genbank_file, 'genbank'):
            # Extract the DEFINITION line from the GenBank record
            definition = record.description
            organism = record.annotations.get('organism', 'Unknown organism')
            print(f"Record definition: {definition}, Organism: {organism}")

            for feature in record.features:
                if feature.type == 'CDS' and 'protein_id' in feature.qualifiers:
                    protein_id = feature.qualifiers['protein_id'][0]
                    
                    # Check if the protein_id is in the mapping
                    if protein_id in accession_map:
                        refseq_accession = accession_map[protein_id]
                        print("Found matching protein ID in GenBank:", protein_id)
                        
                        # Get the start and end of the feature
                        start = feature.location.start
                        end = feature.location.end
                        
                        # Extract upstream and downstream sequences based on strand
                        if feature.location.strand == 1:  # Positive strand
                            upstream_sequence = record.seq[max(0, start - 140):start]
                            downstream_sequence = record.seq[end:end + 531]
                        else:  # Negative strand
                            upstream_sequence = record.seq[end:end + 140].reverse_complement()
                            downstream_sequence = record.seq[max(0, start - 531):start].reverse_complement()

                        # Write to the 5'UTR FASTA file
                        upstream_fasta.write(f">{definition}_upstream\n{upstream_sequence}\n")
                        
                        # Write to the 3'UTR FASTA file
                        downstream_fasta.write(f">{definition}_downstream\n{downstream_sequence}\n")
                        
                        # Collect data for the final Excel file
                        data.append({
                            'Protein Accession': protein_id,
                            'RefSeq Accessions': refseq_accession,
                            'Definition': definition,
                            'Organism': organism,
                            '5\' UTR Sequence': str(upstream_sequence),
                            'CDS Coordinates': f"{start}-{end}",
                            '3\' UTR Sequence': str(downstream_sequence)
                        })
    
    # Write the collected data to an Excel file
    df = pd.DataFrame(data)
    df.to_excel("sequences_and_metadata.xlsx", index=False)
    print("Data has been written to sequences_and_metadata.xlsx")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python second_script.py <genbank_file> <protein_ids_excel_file>")
        sys.exit(1)

    genbank_file = sys.argv[1]
    protein_ids_file = sys.argv[2]
    extract_upstream_downstream_sequence(genbank_file, protein_ids_file)
