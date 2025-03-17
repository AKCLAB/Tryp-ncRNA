import sys
import re
#command: python 3_identify_possible_ncRNA_lncRNA.py transcript_100cov.txt transcript_50cov.txt possible_ncRNA.txt possible_lncRNA.txt
# Inicializar as variáveis globais
i = 1  # numeral do ncRNA
#out_file = "possible_ncRNA.txt" # output file for small ncRNA
#out_file2 = "possible_lncRNA.txt"  # output file for large ncRNA

def process_file(transcript_100, transcript_50, out_file, out_file2):

    with open(transcript_100, 'r') as fh100:
        transcript_lines_100 = fh100.readlines()  # Read transcript_100

    with open(transcript_50, 'r') as fh50:
        transcript_lines_50 = fh50.readlines()  # Read transcript_50
    
        with open(out_file, 'w') as fh_ncRNA, open(out_file2, 'w') as fh_lncRNA:
        # Process the transcripts of transcript_50 (lncRNA)
            for ltra_row in transcript_lines_50:
                ltra_row = ltra_row.strip()  # Remove spaces
                ltra_fields = ltra_row.split("\t")  # Split columns
                ltra_i, ltra_chr, ltra_coordi, ltra_coordf, ltra_strand, ltra_cov = ltra_fields[0], ltra_fields[1], int(ltra_fields[2]), int(ltra_fields[3]), ltra_fields[4], float(ltra_fields[5])
                ltra_length = (ltra_coordf - ltra_coordi) + 1
                # Verify transcript in CDs & not be in the same strand & size
                if ltra_length > 200:
                    fh_lncRNA.write(f"{ltra_chr}\t{ltra_coordi}\t{ltra_coordf}\t{ltra_chr}_l{ltra_i}\t.\t{ltra_strand}\t{ltra_length}\n")

                        # Process the transcripts of transcript_100 (ncRNA)
            for tra_row in transcript_lines_100:
                tra_row = tra_row.strip()  # Remove espace
                tra_fields = tra_row.split("\t")  # split columns
                tra_i, tra_chr, tra_coordi, tra_coordf, tra_strand, tra_cov = tra_fields[0], tra_fields[1], int(tra_fields[2]), int(tra_fields[3]), tra_fields[4], float(tra_fields[5])
                tra_length = (tra_coordf - tra_coordi) + 1
                if 50 <= tra_length <= 200:
                    fh_ncRNA.write(f"{tra_chr}\t{tra_coordi}\t{tra_coordf}\t{tra_chr}_{tra_i}\t.\t{tra_strand}\t{tra_length}\n")


# Import INPUT files
if __name__ == "__main__":
    transcript_100 = sys.argv[1]
    transcript_50 = sys.argv[2]
    out_file = sys.argv[3]
    out_file2 = sys.argv[4]
    process_file(transcript_100, transcript_50, out_file, out_file2)