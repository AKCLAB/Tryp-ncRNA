import pandas as pd
import sys
import numpy as np

#python 9_remake_output.py unique_sort_allncrna.tab df_allncrna.tab df_allncrna.bed df_allncrna.gff

def process_file(sort_allncrna, allncrna_tab, allncrna_bed, allncrna_gff):
    #import all sort aunique ncRNAs
    allncrna = pd.read_csv(sort_allncrna, sep='\t', names=["nc_chr","nc_coordi","nc_coordf","nc_id","nc_strand","nc_length","nc_location","nc_sentido","nc_position"])
    allncrna["nc_clength"] = (allncrna["nc_coordf"] - allncrna["nc_coordi"] + 1) #Calculate the position initial and final of ncRNAa

    ncrna = allncrna[(allncrna["nc_clength"]<=200) & (allncrna["nc_clength"]>=50)] #Select the small ncRNAs
    lncrna = allncrna[allncrna["nc_clength"]>200] #Select the lncRNAs

    chr_ncrna = ncrna["nc_chr"].unique() #List of all chromossomes
    chr_lncrna = lncrna["nc_chr"].unique()

    #Numeration of the ncRNA by chromossome 
    ncrna_list = []
    for uniq_chr in chr_ncrna:
        i = 0
        for allchr in ncrna["nc_chr"]:
            if uniq_chr == allchr:
                i += 1
                ncrna_list.append(allchr+"_ncRNA"+str(i))
    ncrna["ncrna_name"] = ncrna_list
    #Numeration of the lncRNA by chromossome 
    lncrna_list = []
    for uniq_chr in chr_lncrna:
        i = 0
        for allchr in lncrna["nc_chr"]:
            if uniq_chr == allchr:
                i += 1
                lncrna_list.append(allchr+"_lncRNA"+str(i))
    lncrna["ncrna_name"] = lncrna_list

    allncrna_v2 = pd.concat([ncrna, lncrna], axis=0)  # axis=0 means stacking rows
    allncrna_v2.to_csv(allncrna_tab, sep='\t', index=False, header=False)

    #Save format bed of all ncRNAa
    ncrna_bed = allncrna_v2[['nc_chr', 'nc_coordi', 'nc_coordf', 'ncrna_name', 'nc_strand']]
    ncrna_bed.insert(loc=4, column="", value=".")
    ncrna_bed.to_csv(allncrna_bed, sep='\t', index=False, header=False)

    #Create gff format of ncRNA;
    allncrna_v2['nc_clength'] = allncrna_v2['nc_clength'].astype(str)
    ncrna_ncRNA = allncrna_v2[['nc_chr', 'nc_coordi', 'nc_coordf', 'nc_strand']]
    ncrna_ncRNA['info_ncrna'] = "ID="+allncrna_v2['ncrna_name']+":ncRNA;Parent="+allncrna_v2['ncrna_name']+";description=Strand:minus Length:"+allncrna_v2['nc_clength']+" Position:"
    ncrna_ncRNA.insert(loc=1, column ="method", value="InHouseMethodology")
    ncrna_ncRNA.insert(loc=2, column ="molecule", value="ncRNA")
    ncrna_ncRNA.insert(loc=5, column = "", value=".")
    ncrna_ncRNA.insert(loc=7, column = " ", value=".")

    ncrna_exon = allncrna_v2[['nc_chr', 'nc_coordi', 'nc_coordf', 'nc_strand']]
    ncrna_exon['info_ncrna'] = "ID="+allncrna_v2['ncrna_name']+";Parent="+allncrna_v2['ncrna_name']+":ncRNA"
    ncrna_exon.insert(loc=1, column ="method", value="InHouseMethodology")
    ncrna_exon.insert(loc=2, column ="molecule", value="exon")
    ncrna_exon.insert(loc=5, column = "", value=".")
    ncrna_exon.insert(loc=7, column = " ", value=".")

    ncrna_gene = allncrna_v2[['nc_chr', 'nc_coordi', 'nc_coordf', 'nc_strand']]
    ncrna_gene['info_ncrna'] = "ID="+allncrna_v2['ncrna_name']+";description=Strand:minus Length:"+allncrna_v2['nc_clength']+" Position:"
    ncrna_gene.insert(loc=1, column ="method", value="InHouseMethodology")
    ncrna_gene.insert(loc=2, column ="molecule", value="gene")
    ncrna_gene.insert(loc=5, column = "", value=".")
    ncrna_gene.insert(loc=7, column = " ", value=".")
    
    ncrna_peptide = allncrna_v2[['nc_chr', 'nc_coordi', 'nc_coordf', 'nc_strand']]
    ncrna_peptide['info_ncrna'] = "ID="+allncrna_v2['ncrna_name']+";description=Strand:minus Length:"+allncrna_v2['nc_clength']+" Position:"
    ncrna_peptide.insert(loc=1, column ="method", value="InHouseMethodology")
    ncrna_peptide.insert(loc=2, column ="molecule", value="polypeptide")
    ncrna_peptide.insert(loc=5, column = "", value=".")
    ncrna_peptide.insert(loc=7, column = " ", value=".")

    ncrna_gff = pd.concat([ncrna_exon, ncrna_gene, ncrna_ncRNA, ncrna_peptide], axis=0)  # axis=0 means stacking rows
    ncrna_gff = ncrna_gff.sort_values(by=['nc_chr', 'nc_coordi'], ascending=[True, True])
    ncrna_gff.to_csv(allncrna_gff, sep='\t', index=False, header=False)

# Import Files
if __name__ == "__main__":
    sort_allncrna = sys.argv[1]
    allncrna_tab = sys.argv[2]
    allncrna_bed = sys.argv[3]
    allncrna_gff = sys.argv[4]
    process_file(sort_allncrna, allncrna_tab, allncrna_bed, allncrna_gff)