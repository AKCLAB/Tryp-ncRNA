import sys
import argparse

#python3 2_identify_transcript.py "${output_folder}/count_igv.wig" "$output_folder" -threshold 100,50
# Initialize counters for positive and negative strand transcripts
i_pos = 0  # Counter for positive strand transcripts
i_neg = 0  # Counter for negative strand transcripts

def process_file(file, threshold, fh2, out_file):
    global i_pos, i_neg

    with open(file, 'r') as fh:
        chr = None  # Current chromosome
        coordi_pos, coordi_neg = None, None  # Start coordinates for pos and neg
        start_pos, start_neg = 0, 0  # Flags to indicate open ranges
        prev_coord = None

        for line_number, row in enumerate(fh):
            # Skip the first two lines
            if line_number < 2:
                continue
            if row.startswith("variableStep"):
                # Handle chromosome transition: close open ranges
                if start_pos == 1:
                    fh2.write(f"ncRNA{i_pos}\t{chr}\t{coordi_pos}\t{prev_coord}\t+\t{threshold}\n")
                    start_pos = 0
                if start_neg == 1:
                    fh2.write(f"ncRNA{i_neg}\t{chr}\t{coordi_neg}\t{prev_coord}\t-\t{threshold}\n")
                    start_neg = 0

                # Update chromosome name
                row = row.strip()
                name = row.split(" ")
                chr = name[1].split("=")[1]
                continue

            # Process data rows
            fields = row.strip().split("\t")
            if len(fields) < 3:  # Ensure there are enough columns
                continue
            coord = fields[0]
            qtd_pos = float(fields[1])  # Column for "pos"
            qtd_neg = float(fields[2])  # Column for "neg"
            # Process positive strand
            if qtd_pos >= threshold and start_pos == 0:
                start_pos = 1
                coordi_pos = coord  # Save start coordinate

                i_pos += 1
            elif start_pos == 1 and qtd_pos < threshold:
                fh2.write(f"ncRNA{i_pos}\t{chr}\t{coordi_pos}\t{prev_coord}\t+\t{threshold}\n")
                start_pos = 0
                coordi_pos = None
                #IF coordi_pos = prev_coord + OR +1 OR +2 
                #fh2.write(f"ncRNA{i_pos}\t{chr}\t{coordi_pos}\t{prev_coord}\t+\t{threshold}\n")

            # Process negative strand
            if qtd_neg >= threshold and start_neg == 0:
                start_neg = 1
                coordi_neg = coord  # Save start coordinate
                i_neg += 1
            elif start_neg == 1 and qtd_neg < threshold:
                fh2.write(f"ncRNA{i_neg}\t{chr}\t{coordi_neg}\t{prev_coord}\t-\t{threshold}\n")
                start_neg = 0
                coordi_neg = None

            prev_coord = coord

        # Handle open ranges at the end of the file
        if start_pos == 1:
            fh2.write(f"ncRNA{i_pos}\t{chr}\t{coordi_pos}\t{prev_coord}\t+\t{threshold}\n")
        if start_neg == 1:
            fh2.write(f"ncRNA{i_neg}\t{chr}\t{coordi_neg}\t{prev_coord}\t-\t{threshold}\n")

if __name__ == "__main__":
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Process a transcript file with thresholds.")
    parser.add_argument("file", type=str, help="Input file (e.g., count_igv.wig)")
    parser.add_argument("output_directory", type=str, help="Output directory")
    parser.add_argument("-threshold", type=str, required=True, help="Comma-separated list of thresholds (e.g., 100,50)")

    # Parse the arguments
    args = parser.parse_args()

    # Parse the threshold argument into a list of integers
    thresholds = [int(x) for x in args.threshold.split(",")]

    # Ensure the output directory ends with a slash
    output_directory = args.output_directory
    if not output_directory.endswith("/"):
        output_directory += "/"

    # Process the input file for each threshold
    for threshold in thresholds:
        # Create an output file for each threshold
        out_file = f"{output_directory}transcript_{threshold}cov.txt"
        with open(out_file, 'w') as fh2:
            # Process the file with the given threshold
            process_file(args.file, threshold, fh2, out_file)