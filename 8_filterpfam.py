import pandas as pd
import sys
from Bio import SeqIO

#python3 8_filterpfam.py ncrna_pfam-cov80-max1.tab header_to_remove.tsv all_ncrna.fasta filtered_pfam.fasta
#Name files
#blast_tab = "ncrna_pfam-cov80-max1.tab"
#header_remove = header_to_remove.tsv
#all_ncrna = "all_ncrna.fasta"
#filtered_pfam = "filtered_pfam.fasta"

def process_file(blast_tab,header_remove, all_ncrna, filtered_pfam):

    #import blast output with header names
    blast_pfam = pd.read_csv(blast_tab, sep='\t', names=["query id","subject id","identity","query coverage per subject","alignment length","query length","subject length","q. start","q. end","s. start","s. end","evalue","bit score", "subject title"])
    blast_pfam2 = blast_pfam[blast_pfam["identity"] > 90] #filter with identity > 90%
    remove = blast_pfam2["query id"]
    remove.to_csv(header_remove, sep='\t', index=False, header=False) #save the ncRNA blast with pfam database

    headers_to_remove = set(line.strip() for line in remove) #Go ID  by ID

    # Filter the fasta sequences that not coding region
    with open(all_ncrna) as infile, open(filtered_pfam, "w") as outfile:
        for record in SeqIO.parse(infile, "fasta"):
            if record.id not in headers_to_remove:
                SeqIO.write(record, outfile, "fasta")

# Import Files
if __name__ == "__main__":
    blast_tab = sys.argv[1]
    header_remove = sys.argv[2]
    all_ncrna = sys.argv[3]
    filtered_pfam = sys.argv[4]
    process_file(blast_tab, header_remove, all_ncrna, filtered_pfam)