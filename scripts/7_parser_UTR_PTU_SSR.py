import sys
import re
#command: python 7_parser_UTR_sense_antisense.py annotation_ncRNA_final.txt annotation_lncRNA.txt possible_ptu.txt possible_ssr.txt ncRNAs_location_direction.bed LBRAZ_M2903.Dec2022.gff UTRme_fiveutr.tsv UTRme_threeutr.tsv
#out_file = "ncRNAs_location_direction.bed" # output file with location, position, sense of noncoding RNA
intergenic = "intergenic"
# Define the function to extract the ID value from the 9th column
def extract_info(text):
    id_match = re.search(r'ID=([^;]+)', text)
    id_value = id_match.group(1) if id_match else None
    return id_value
#
def process_file(ncrna, lcrna, ptu, ssr, out_file, gff, sl=None, polya=None):
    #global out_file
    with open(ncrna, 'r') as file1, open(lcrna, 'r') as file2,  open(gff, 'r') as fh_cds: #open ncRNA and lcRNA
        ncrna_lines = file1.readlines() #read all lines
        lcrna_lines = file2.readlines()
        all_rna_lines = ncrna_lines + lcrna_lines  # Join lines of two files  
        cds_lines = []
        for cds_row in fh_cds:
            if cds_row.startswith("#"):
                continue
            # Split the line into columns
            fields = cds_row.strip().split("\t")
            # Check if the line is a CDS entry
            if fields[2] in ["CDS", "exon"]:  
                cds_lines.append(fields)


    with open(ptu, 'r') as file3, open(ssr, 'r') as file6:
        #ptus_lines = file3.readlines() #reads from policistronic coordenates lines
        sense_lines = file3.readlines() + file6.readlines() #reads from policistronic coordenates lines
    utr5_lines = []
    #utr3_lines = []

    if sl:
        with open(sl, 'r') as file4: # read leader sequence and polya position
            utr5_lines = file4.readlines()
    if polya:
        with open(polya, 'r') as file5:
            utr3_lines = file5.readlines()

    with open(out_file, 'w') as ncRNA: #write output file
            
        #Read transcript file
        for row in all_rna_lines:
            row = row.strip()  # Remove spaces
            fields = row.split("\t")  # Split columns
            nc_chr, nc_coordi, nc_coordf, nc_id, nc_strand, nc_length = fields[0], int(fields[1]), int(fields[2]), fields[3], fields[5], fields[6]
            nc_sentido = "" #empty list to save the direction of PTU and ncRNA
            # Verify ncRNA in PTUs
            for sense_row in sense_lines:
                sense_row = sense_row.strip() # Remove spaces
                sense_fields = sense_row.split("\t") # Split columns
                sense_chr, sense_strand, sense_coordi, sense_coordf, sense_name= sense_fields[0], sense_fields[1], int(sense_fields[2]), int(sense_fields[3]), sense_fields[4]
                if nc_chr == sense_chr: # If be in the same chromossome and overlapping ptu regions
                    if sense_fields[4] == "PTU":
                        if((nc_coordi >= sense_coordi and nc_coordi< sense_coordf) or 
                        (nc_coordf > sense_coordi and nc_coordf <= sense_coordf) or 
                        (nc_coordi >= sense_coordi and nc_coordf <= sense_coordf)):
                            if nc_strand == sense_strand: #if the ncRNA be in the same strand
                                nc_sentido = "sense"
                            else: #if the ncRNA NOT TO be in the same strand
                                nc_sentido = "antisense" 
                            break
                    else:
                        if((nc_coordi >= sense_coordi and nc_coordi< sense_coordf) or 
                        (nc_coordf > sense_coordi and nc_coordf <= sense_coordf) or 
                        (nc_coordi >= sense_coordi and nc_coordf <= sense_coordf)):
                            nc_sentido = sense_name

                            break

            nc_position = "" 
            matched = False
            for cds_fields in cds_lines:
                cds_id = extract_info(cds_fields[8])
                #print(cds_id)
                cds_chr, cds_coordi, cds_coordf, cds_strand = cds_fields[0],  int(cds_fields[3]), int(cds_fields[4]), cds_fields[6]


                # Verify transcript in CDs & not be in the same strand & size
                if (nc_chr == cds_chr and 
                    ((nc_coordi >= cds_coordi and nc_coordi <= cds_coordf) or (nc_coordf >= cds_coordi and nc_coordf <= cds_coordf) or
                    (cds_coordi >= nc_coordi and cds_coordi <= nc_coordf) or
                    (cds_coordf >= nc_coordi and cds_coordf <= nc_coordf))):

                    if nc_strand != cds_strand:
                        #prin\t(tra_chr, tra_coordi, tra_coordf, tra_strand, tra_length, cds_id)
                        nc_position = cds_id
                        #matched = True
                        #break  # Match found for each loop
                    else:  # Mesmo strand, não é intergênico
                        nc_position = intergenic
                    matched = True
                    break
            if not matched:        
                for utr5_line in utr5_lines:
                    # Split the line into columns
                    utr5_fields = utr5_line.strip().split("\t")
                    # Check if the line is a CDS entry
                    # Extract the ID from the 9th column (index 8)
                    utr5_chr, utr5_coordi, utr5_coordf, utr5_strand, utr5_name = utr5_fields[0], int(utr5_fields[1]), int(utr5_fields[2]), utr5_fields[3], utr5_fields[5]
                    #print(utr5_name, utr5_coordi)
                    if nc_chr == utr5_chr and nc_strand == utr5_strand: # If be in the same chromossome and overlapping UTR regions
                        if((nc_coordi > utr5_coordi and nc_coordi < utr5_coordf) or (nc_coordf > utr5_coordi and nc_coordf < utr5_coordf) or (nc_coordi< utr5_coordi and nc_coordf > utr5_coordf)):
                            nc_position = utr5_name
                            matched = True
                            break
            if not matched:
                for utr3_line in utr3_lines:
                    utr3_fields = utr3_line.strip().split("\t")
                    utr3_chr, utr3_coordi, utr3_coordf, utr3_strand, utr3_name = utr3_fields[0], int(utr3_fields[1]), int(utr3_fields[2]), utr3_fields[3], utr3_fields[5]
                    if nc_chr == utr3_chr and nc_strand == utr3_strand: # If be in the same chromossome and overlapping UTR regions
                        if((nc_coordi > utr3_coordi and nc_coordi < utr3_coordf) or (nc_coordf > utr3_coordi and nc_coordf < utr3_coordf) or (nc_coordi< utr3_coordi and nc_coordf > utr3_coordf)):
                            nc_position = utr3_name
                            matched = True
                            break

            if not matched:
                nc_position = intergenic

            output_line = (
                f"{nc_chr}\t{nc_coordi}\t{nc_coordf}\t{nc_id}\t.\t{nc_strand}\t{nc_length}\t" # SAVE THE COLUNS 
                f"{nc_sentido if nc_sentido else 'undetermined'}\t{nc_position if nc_position else 'undetermined'}\n"
            ) #ADD COLUNS WITH INFORMATION
            ncRNA.write(output_line)       
# Import INPUT files
if __name__ == "__main__":
    ncrna = sys.argv[1]
    lcrna = sys.argv[2]
    ptu = sys.argv[3]
    ssr = sys.argv[4]
    out_file = sys.argv[5]
    gff = sys.argv[6]
    sl = sys.argv[7] if len(sys.argv) > 7 else None
    polya = sys.argv[8] if len(sys.argv) > 8 else None
    process_file(ncrna, lcrna, ptu, ssr, out_file, gff, sl, polya)
    #process_file(ncrna, lcrna, out_file, gff, sl, polya)