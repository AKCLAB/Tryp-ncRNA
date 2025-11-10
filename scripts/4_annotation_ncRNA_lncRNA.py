import sys
import re
#command: python 4_annotation_ncRNA_lncRNA.py possible_ncRNA.txt possible_lncRNA.txt LBRAZ_M2903.Dec2022.gff annotation_ncRNA.txt annotation_lncRNA.txt

# Inicializar as variáveis globais
i = 1  # numeral do ncRNA
#out_file = "annotation_ncRNA.txt" # output file for small ncRNA
#out_file2 = "annotation_lncRNA.txt"  # output file for large ncRNA
intergenic = "intergenic"  # name for intergenic position

# Function to extract the ID from the 9th column in a GFF file
def extract_info(text):
    id_match = re.search(r'ID=([^;]+)', text)  # extract ID after "ID=" and before ";"
    id_value = id_match.group(1) if id_match else None
    return id_value

def process_file(transcript_ncrna, transcript_lncrna, gff, out_file, out_file2):
    global i, intergenic

    with open(transcript_ncrna, 'r') as fh_ncRNA:
        transcript_ncrna = fh_ncRNA.readlines()  # Read transcript_100cov

    with open(transcript_lncrna, 'r') as fh_lncRNA:
        transcript_lncrna = fh_lncRNA.readlines()  # Read transcript_50cov

    with open(gff, 'r') as fh_cds:
        # Write in output files
        with open(out_file, 'w') as fh_ncRNA, open(out_file2, 'w') as fh_lncRNA:
            # Process the large transcripts of transcript_50cov (lncRNA)
            for ltra_row in transcript_lncrna:
                fh_cds.seek(0) 
                ltra_row = ltra_row.strip()  # Remove spaces
                ltra_fields = ltra_row.split("\t")  # Split columns
                ltra_chr, ltra_coordi, ltra_coordf, ltra_i, ltra_strand, ltra_length = ltra_fields[0],  int(ltra_fields[1]), int(ltra_fields[2]), ltra_fields[3], ltra_fields[5], int(ltra_fields[6])
                overlap_same_strand = False
                overlap_diff_strand = False
                # Verify lncRNA overlapping on CDSs
                for cds_row in fh_cds:
                    if cds_row.startswith("#"):
                        continue
                    # Split the line into columns
                    cds_fields = cds_row.strip().split("\t")
                    # Check if the line is a CDS entry
                    if cds_fields[2] in ["CDS", "exon"]:
                        #cds_id = extract_info(cds_fields[8])
                        cds_chr, cds_coordi, cds_coordf, cds_strand = cds_fields[0],  int(cds_fields[3]), int(cds_fields[4]), cds_fields[6]
                        # Verify if the transcript is overlapping CDSs with the same strand location
                        if (ltra_chr == cds_chr and 
                            ((ltra_coordi >= cds_coordi and ltra_coordi <= cds_coordf) or (ltra_coordf >= cds_coordi and ltra_coordf <= cds_coordf) or
                            (cds_coordi >= ltra_coordi and cds_coordi <= ltra_coordf) or
                            (cds_coordf >= ltra_coordi and cds_coordf <= ltra_coordf))):
                            # Mark if it's the same or different strand
                            if ltra_strand != cds_strand:
                                overlap_diff_strand = True
                            else:
                                overlap_same_strand = True
                # Keep it (only if all overlays are with different strands)
                if overlap_diff_strand and not overlap_same_strand:
                    fh_lncRNA.write(f"{ltra_chr}\t{ltra_coordi}\t{ltra_coordf}\t{ltra_i}\t.\t{ltra_strand}\t{ltra_length}\n")
                elif not overlap_diff_strand and not overlap_same_strand:
                # It does not overlap with any CDS → it is an intergenic lncRNA.
                    fh_lncRNA.write(f"{ltra_chr}\t{ltra_coordi}\t{ltra_coordf}\t{ltra_i}\t.\t{ltra_strand}\t{ltra_length}\n")

            # Process the transcripts of transcript_100 (ncRNA)
            for tra_row in transcript_ncrna:
                fh_cds.seek(0) 
                tra_row = tra_row.strip()  # Remove espace
                tra_fields = tra_row.split("\t")  # split columns
                tra_chr, tra_coordi, tra_coordf, tra_i, tra_strand, tra_length = tra_fields[0], int(tra_fields[1]), int(tra_fields[2]), tra_fields[3], tra_fields[5], int(tra_fields[6])
                matched = False  # Control variable to verify match
                overlap_same_strand = False
                overlap_diff_strand = False
                # Verify lncRNA overlapping on CDsA
                for cds_row in fh_cds:
                    if cds_row.startswith("#"):
                        continue
                    # Split the line into columns
                    cds_fields = cds_row.strip().split("\t")
                    # Check if the line is a CDS entry
                    if cds_fields[2] in ["CDS", "exon"]:
                        cds_chr, cds_coordi, cds_coordf, cds_strand = cds_fields[0],  int(cds_fields[3]), int(cds_fields[4]), cds_fields[6]
                        # Verify transcript in CDs & not be in the same strand & size
                        if tra_chr == cds_chr and (((tra_coordi >= cds_coordi and tra_coordi <= cds_coordf) or
                            (tra_coordf >= cds_coordi and tra_coordf <= cds_coordf) or
                            (cds_coordi >= tra_coordi and cds_coordi <= tra_coordf) or
                            (cds_coordf >= tra_coordi and cds_coordf <= tra_coordf))):
                            # Mark if it's the same or different strand
                            if tra_strand != cds_strand:
                                overlap_diff_strand = True
                            else: 
                                overlap_same_strand = True

                # Keep it (only if all overlays are with different strands)
                if overlap_diff_strand and not overlap_same_strand:
                    fh_ncRNA.write(f"{tra_chr}\t{tra_coordi}\t{tra_coordf}\t{tra_i}\t.\t{tra_strand}\t{tra_length}\n")
                elif not overlap_diff_strand and not overlap_same_strand:
                # It does not overlap with any CDS → it is an intergenic lncRNA.
                    fh_ncRNA.write(f"{tra_chr}\t{tra_coordi}\t{tra_coordf}\t{tra_i}\t.\t{tra_strand}\t{tra_length}\n")

# Import INPUT files
if __name__ == "__main__":
    transcript_ncrna = sys.argv[1]
    transcript_lncrna = sys.argv[2]
    gff = sys.argv[3]
    out_file = sys.argv[4]
    out_file2 = sys.argv[5]
    process_file(transcript_ncrna, transcript_lncrna, gff, out_file, out_file2)
