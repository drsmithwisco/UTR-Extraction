# UTR-Extraction
# Workflow for extraction of the 5' and 3' untranslated regions of a protein sequence 

Dependencies to download: 
- Use pip to install these dependencies if they are not installed  already. 

pip install pandas biopython openpyxl

- clone the repository into a directory of your choice: 

git clone https://github.com/drsmithwisco/UTR-Extraction

- navigate to the directory and open up the scripts subdirectory. 
This contains three python scripts. You will also need to have the excel file 
containing the protein accessions with the column header exactly called Protein Accession

-You must have the excel file containing the protein refseqs in the same directory as the python scripts

- Make sure to open up the ExtractRefseq_from_IPG.py and the Refseq_to_Genbanks.py 
python scripts and add your email address for NCBI. 

Usage:
python3 ExtractRefseq_from_IPG.py

General Workflow:
The first script reads an excel file and searches the protein accession IDs against the Identical protein groups database 
Then it extracts the refseq associated with it. 

The second script extracts the genbank file per refseq accession and concatenates all of the genbanks into a single file. 

The third script locates the protein accessions in the genbank and extracts nucleotide
upstream and downstream of this and creates a final excel file with the metadata,
