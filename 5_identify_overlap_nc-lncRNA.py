import sys
#command: python 5_identify_overlap_nc-lncRNA.py annotation_ncRNA.txt annotation_lncRNA.txt annotation_ncRNA_final.txt
# Inicializar as variáveis globais
#out_file = "annotation_ncRNA_final.txt"

def process_file(possible_ncRNA, possible_lncRNA, out_file):
    #global out_file

    with open(possible_ncRNA, 'r') as ncRNA:
        transcript_nc = ncRNA.readlines()  # Read ncRNA file

    with open(possible_lncRNA, 'r') as lncRNA:
        transcript_lnc = lncRNA.readlines()  # Read lncRNA file

    # Open output file for writing filtered results
    with open(out_file, 'w') as fh_ncRNA:
        for nc_row in transcript_nc:
            nc_row = nc_row.strip()
            nc_fields = nc_row.split("\t")
            nc_chr, nc_coordi, nc_coordf, nc_strand = nc_fields[0], int(nc_fields[1]), int(nc_fields[2]), nc_fields[5]

            keep_ncRNA = True  # Default to keeping ncRNA unless overlap invalidates it

            for lnc_row in transcript_lnc:
                lnc_row = lnc_row.strip()
                lnc_fields = lnc_row.split("\t")
                lnc_chr, lnc_coordi, lnc_coordf, lnc_strand = lnc_fields[0], int(lnc_fields[1]), int(lnc_fields[2]), lnc_fields[5]

                if nc_chr == lnc_chr:  # Check chromosome match
                    # Check for overlap
                    if (nc_coordi <= lnc_coordf and nc_coordf >= lnc_coordi):
                        # If overlap exists and the strands match, mark for exclusion
                        if nc_strand == lnc_strand:
                            keep_ncRNA = False
                            break  # No need to check further

            # Write ncRNA to file if it should be kept
            if keep_ncRNA:
                fh_ncRNA.write(nc_row + "\n")

# Import input files
if __name__ == "__main__":
    possible_ncRNA = sys.argv[1]
    possible_lncRNA = sys.argv[2]
    out_file = sys.argv[3]
    process_file(possible_ncRNA, possible_lncRNA, out_file)
