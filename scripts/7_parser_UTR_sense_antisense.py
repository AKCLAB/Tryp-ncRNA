import sys
#command: python 7_parser_UTR_sense_antisense.py annotation_ncRNA_final.txt annotation_lncRNA.txt possible_ptu.txt ncRNAs_location_direction.bed UTRme_fiveutr.tsv UTRme_threeutr.tsv
#out_file = "ncRNAs_location_direction.bed" # output file with location, position, sense of noncoding RNA

def process_file(ncrna, lcrna, ptu, out_file, sl=None, polya=None):
    #global out_file
    with open(ncrna, 'r') as file1, open(lcrna, 'r') as file2: #open ncRNA and lcRNA
        ncrna_lines = file1.readlines() #read all lines
        lcrna_lines = file2.readlines()
        all_rna_lines = ncrna_lines + lcrna_lines  # Join lines of two files    
        with open(ptu, 'r') as file3:
            ptus_lines = file3.readlines() #reads from policistronic coordenates lines
        utr5_lines = []
        utr3_lines = []
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
                        nc_chr, nc_coordi, nc_coordf, nc_id, nc_strand, nc_length, nc_location = fields[0], int(fields[1]), int(fields[2]), fields[3], fields[5], fields[6], fields[7]
                        nc_sentido = "" #empty list to save the direction of PTU and ncRNA
                        nc_position = "" #empty list to save the position of ncRNA in 5UTR or 3UTR
                        # Verify ncRNA in PTUs
                        for ptus_row in ptus_lines:
                            ptus_row = ptus_row.strip() # Remove spaces
                            ptus_fields = ptus_row.split("\t") # Split columns
                            ptus_chr, ptus_coordi, ptus_coordf, ptus_strand = ptus_fields[0], int(ptus_fields[2]), int(ptus_fields[3]), ptus_fields[1]
                            if nc_chr == ptus_chr: # If be in the same chromossome and overlapping ptu regions
                                if((nc_coordi >= ptus_coordi and nc_coordi< ptus_coordf) or 
                                (nc_coordf > ptus_coordi and nc_coordf <= ptus_coordf) or 
                                (nc_coordi >= ptus_coordi and nc_coordf <= ptus_coordf)):
                                    if nc_strand == ptus_strand: #if the ncRNA be in the same strand
                                        nc_sentido = "sense"
                                    else: #if the ncRNA NOT TO be in the same strand
                                        nc_sentido = "antisense" 
                                    break
                        # Verify lncRNA overlapping on CDsA
                        for utr5_line in utr5_lines:
                            # Split the line into columns
                            utr5_fields = utr5_line.strip().split("\t")
                            # Check if the line is a CDS entry
                            # Extract the ID from the 9th column (index 8)
                            utr5_chr, utr5_coordi, utr5_coordf, utr5_name = utr5_fields[0], int(utr5_fields[1]), int(utr5_fields[2]), utr5_fields[5]
                            #print(utr5_name, utr5_coordi)
                            if nc_chr == utr5_chr: # If be in the same chromossome and overlapping UTR regions
                                if((nc_coordi > utr5_coordi and nc_coordi < utr5_coordf) or 
                                    (nc_coordf > utr5_coordi and nc_coordf < utr5_coordf) or 
                                    (nc_coordi< utr5_coordi and nc_coordf > utr5_coordf)):
                                    nc_position = utr5_name
                                    break
                        for utr3_line in utr3_lines:
                            utr3_fields = utr3_line.strip().split("\t")
                            #utr3_id = extract_info(utr3_fields[8])
                            utr3_chr, utr3_coordi, utr3_coordf, utr3_name = utr3_fields[0], int(utr3_fields[1]), int(utr3_fields[2]), utr3_fields[5]
                            if nc_chr == utr3_chr: # If be in the same chromossome and overlapping UTR regions
                                if((nc_coordi > utr3_coordi and nc_coordi < utr3_coordf) or (nc_coordf > utr3_coordi and nc_coordf < utr3_coordf) or (nc_coordi< utr3_coordi and nc_coordf > utr3_coordf)):
                                    nc_position = utr3_name
                                    break
                        output_line = (
                            f"{nc_chr}\t{nc_coordi}\t{nc_coordf}\t{nc_id}\t.\t{nc_strand}\t{nc_length}\t{nc_location}\t" # SAVE THE COLUNS 
                            f"{nc_sentido if nc_sentido else 'undetermined'}\t{nc_position if nc_position else 'undetermined'}\n"
                        ) #ADD COLUNS WITH INFORMATION
                        ncRNA.write(output_line)       
# Import INPUT files
if __name__ == "__main__":
    ncrna = sys.argv[1]
    lcrna = sys.argv[2]
    ptu = sys.argv[3]
    out_file = sys.argv[4]
    sl = sys.argv[5] if len(sys.argv) > 4 else None
    polya = sys.argv[6] if len(sys.argv) > 5 else None
    
    process_file(ncrna, lcrna, ptu, out_file, sl, polya)

