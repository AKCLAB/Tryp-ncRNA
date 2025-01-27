import pandas as pd
import sys

#code: python3 8_select_ncrna.py ncRNAs_location_direction.bed ncrna_pfam-cov80-max1.tab final_ncRNAs_location_direction.tab

#input & output files
#ncrna_bed = "ncRNAs_location_direction.bed"
#blast_tab = ncrna_pfam-cov80-max1.tab
#output= "final_ncRNAs_location_direction.bed"
# Define the function

def process_file(ncrna_file, blast_tab, final_ncrna):

    # Importar o arquivo usando pandas
    ncrna = pd.read_csv(ncrna_file, sep='\t', names=["nc_chr", "nc_coordi", "nc_coordf", "nc_id", "", "nc_strand","nc_length", "nc_location", "nc_sentido", "nc_position"], dtype={"nc_coordi": "str", "nc_coordf": "str"})

    # Combine the columns with separators
    ncrna['info_ncrna'] = ncrna['nc_id'] + "::" + ncrna['nc_chr'] + ":" + ncrna['nc_coordi'] + "-" + ncrna['nc_coordf']

    #column_removed = pd.read_csv(removed_file, names=["names"]) #list of ncRNA with match in pfam proteins
    blast_pfam = pd.read_csv(blast_tab, sep='\t', names=["query id","subject id","identity","query coverage per subject","alignment length","query length","subject length","q. start","q. end","s. start","s. end","evalue","bit score", "subject title"])
    
    # From the list all predicted ncRNA remove the included in removed_file
    filtered_df = ncrna[~ncrna['info_ncrna'].isin(blast_pfam["query id"])]
    filtered_df.drop('info_ncrna', axis=1, inplace=True) 
    filtered_df.to_csv(final_ncrna, sep='\t', index=False, header=False) # final predicted ncRNA

# Import Files
if __name__ == "__main__":
    ncrna_file = sys.argv[1]
    blast_tab = sys.argv[2]
    final_ncrna = sys.argv[3]
    process_file(ncrna_file, blast_tab, final_ncrna)
