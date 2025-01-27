#!/bin/bash
# Function to display usage information
usage() {
    echo "Usage: $0 -dir_list <file> -output <path> -db <path> -threads <number> -reffasta <file> -refgff <file> -utr5 <file> -utr3 <file>"
    echo "  -dir_list: File containing a list of directories to process"
    echo "  -output: Base directory where outputs will be stored"
    echo "  -db: Base directory to dave Pfam database"
    echo "  -threads: Number of threads to use"
    echo "  -reffasta: Reference genome FASTA file"
    echo "  -refgff: Annotation GFF file"
    echo "  -utr5: Annotation file for 5' UTRs (optional)"
    echo "  -utr3: Annotation file for 3' UTRs (optional)"
    exit 1
}

#bash pipe_ncrna.sh -dir_list listdir.txt -output /path/to/TriTry-ncRNA -db /path/to/TriTry-ncRNA/db -threads 20 -reffasta /path/to/work/TriTrypDB-30_LbraziliensisMHOMBR75M2903_Genome.fasta -refgff /path/to/work/TriTrypDB-30_LbraziliensisMHOMBR75M2903.gff -utr5 /path/to/work/UTRme_fiveutr.tsv -utr3 /path/to/work/UTRme_threeutr.tsv
#bash pipe_ncrna.sh -dir_list /path/to/TriTry-ncRNA/pro -output /path/to/TriTry-ncRNA/pro_out -db -threads 20 -fasta /path/to/work/TriTrypDB-30_LbraziliensisMHOMBR75M2903_Genome.fasta -gff /path/to/work/TriTrypDB-30_LbraziliensisMHOMBR75M2903.gff  -utr5 /path/to/work/UTRme_fiveutr.tsv -utr3 /path/to/work/UTRme_threeutr.tsv

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        -dir_list)
        dir_list="$2"
        shift # past argument
        shift # past value
        ;;
        -output)
        output_base="$2"
        shift
        shift
        ;;
        -db)
        database="$2"
        shift
        shift
        ;;
        -threads)
        threads="$2"
        shift
        shift
        ;;
        -reffasta)
        fasta="$2"
        shift
        shift
        ;;
        -refgff)
        gff="$2"
        shift
        shift
        ;;
        -utr5)
        utr5="$2"
        shift
        shift
        ;;
        -utr3)
        utr3="$2"
        shift
        shift
        ;;
        *)    # unknown option
        usage
        ;;
    esac
done

# Check if all variables are set
if [ -z "$dir_list" ] || [ -z "$output_base" ] || [ -z "$database" ] || [ -z "$threads" ] || [ -z "$fasta" ]; then
    usage
fi

# Now use these variables in your script
echo "FASTQ Directory: $dir_list"
echo "Output Folder: $output_base"
echo "Pfam database: $database"
echo "Threads: $threads"
echo "Reference FASTA: $fasta"
echo "GFF File: $gff"
echo "optional UTR5 File: $utr5"
echo "optional UTR3 File: $utr3"

#verify the index files or run reference index 
ref_name=$(basename "$fasta" .fasta)
# Verificar si no existen archivos .bt2
if ! ls "$ref_name"*.bt2 1> /dev/null 2>&1; then
    echo "Running bowtie-build"
    bowtie2-build "${ref_name}.fasta" "$ref_name"
    echo "Index created successfully"
else
    echo "BT2 files found, skipping bowtie-build"
fi

# Verify the index file
if ! ls "$ref_name"*.fai 1> /dev/null 2>&1; then
    echo "Running bowtie-build"
    samtools faidx "${ref_name}.fasta"
    echo "Index created successfully"
else
    echo "BT2 files found, skipping samtools"
fi

echo "Verify the Pfam database"
if [ ! -f "$database/pfam_database.dmnd" ]; then
    mkdir -p "$database"
    wget -O "$database/Pfam-A.fasta.gz" https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.fasta.gz
    gunzip -c "$database/Pfam-A.fasta.gz" > "$database/Pfam-A.fasta"
    makeblastdb -in "$database/Pfam-A.fasta" -dbtype prot -out "$database/pfam_database"
    rm "$database/Pfam-A.fasta.gz"
fi
echo "Pfam database created successfully"

# Process each directory listed in the input file
while IFS= read -r subdir; do
    # Skip empty lines
    [[ -z "$subdir" ]] && continue

    # Check if the input directory exists
    if [ ! -d "$subdir" ]; then
        echo "Warning: Directory '$subdir' not found. Skipping."
        continue
    fi

    # Extract the base name of the input directory
    dir_name=$(basename "$subdir")

    # Construct the output directory name
    output_folder="${output_base}/${dir_name}_out"
    # Ensure the output directory exists
    #mkdir -p "$output_folder"
    if [ ! -d "$output_folder" ]; then
        mkdir "$output_folder"
    fi

    name_samples=()
    # Loop por todas as amostras
    for fastq_file in "${subdir}"/*_1.fastq.gz; do
        #Nome do arquivo .fastq (sem a extensão)
        fastq_name=$(basename "$fastq_file" _1.fastq.gz)
        name_samples+=("$fastq_name")
    done
    echo "${name_samples}"

    # Run fastqc
    #fastqc * --threads "$threads"
    echo "File fastqc successfully"

    echo "Running Bowtie2"
    # Loop through each sample

    for sample in "${name_samples[@]}"; do
        echo "Processing $sample..."

        # Check if the sorted BAM and index files already exist
        if [ -f "${output_folder}/mapped_${sample}_sorted.bam" ] && [ -f "${output_folder}/mapped_${sample}_sorted.bai" ]; then
            echo "Files for $sample already exist. Skipping Bowtie2 for this sample."
            continue
        fi
        # Run bowtie2
        bowtie2 \
        -N 1 \
        -p "$threads" \
        --local \
        -x "$ref_name" \
        -1 "${subdir}/${sample}_1.fastq.gz" \
        -2 "${subdir}/${sample}_2.fastq.gz" \
        -S "${output_folder}/mapped_${sample}.sam" 2> "${output_folder}/mapped_${sample}.log"

        # Convert SAM to BAM using samtools
        samtools view -bS "${output_folder}/mapped_${sample}.sam" > "${output_folder}/mapped_${sample}.bam"
        rm "${output_folder}/mapped_${sample}.sam"
        samtools sort "${output_folder}/mapped_${sample}.bam" -o "${output_folder}/mapped_${sample}_sorted.bam"
        rm "${output_folder}/mapped_${sample}.bam"
        samtools index "${output_folder}/mapped_${sample}_sorted.bam" "${output_folder}/mapped_${sample}_sorted.bai"
        echo "Mapping done successfully for ${sample}"
    done
    echo "Processing complete."

    # String withs name sample separated by spaces and -I
    inputs=""
    for sample in "${name_samples[@]}"; do
        inputs+=" I=${output_folder}/mapped_${sample}_sorted.bam"
    done

    echo "Running Picard"
    picard MergeSamFiles ${inputs} USE_THREADING=true O="${output_folder}/transcript_all.sorted.merged_files.bam"
    picard BuildBamIndex I="${output_folder}/transcript_all.sorted.merged_files.bam"
    #echo "File picard successfully"

    echo "Running igvtools"
    igvtools count --strands second --windowSize 1 "${output_folder}/transcript_all.sorted.merged_files.bam" "${output_folder}/count_igv.wig,count_igv.tdf" "$fasta"
    echo "File igvtools successfully"

    echo "Running 2_identify_transcript.py"
    # Call 2_identify_transcript.py script
    python3 2_identify_transcript.py "${output_folder}/count_igv.wig" "$output_folder" -threshold 100,50
    echo "Identify transcript done successfully"

    echo "Running 3_identify_possible_ncRNA_lncRNA.py"
    # Call 3_identify_possible_ncRNA_lncRNA.py script
    python3 3_identify_possible_ncRNA_lncRNA.py "${output_folder}/transcript_100cov.txt" "${output_folder}/transcript_50cov.txt" "${output_folder}/possible_ncRNA.txt" "${output_folder}/possible_lncRNA.txt"
    echo "Identify ncRNA and lncRNA done successfully"

    echo "Running 4_annotation_ncRNA_lncRNA.py script"
    # Call 4_annotation_ncRNA_lncRNA.py script
    python3 4_annotation_ncRNA_lncRNA.py "${output_folder}/possible_ncRNA.txt" "${output_folder}/possible_lncRNA.txt" "$gff" "${output_folder}/annotation_ncRNA.txt" "${output_folder}/annotation_lncRNA.txt"
    #merge of cycle change
    echo "Position annnotation of ncRNA and lncRNA"

    echo "Running 5_identify_overlap_nc-lncRNA.py script"
    # Call 5_identify_overlap_nc-lncRNA.py script
    python3 5_identify_overlap_nc-lncRNA.py "${output_folder}/annotation_ncRNA.txt" "${output_folder}/annotation_lncRNA.txt" "${output_folder}/annotation_ncRNA_final.txt"
    echo "Identify ncRNA and lncRNA do not overlapping"

    echo "Running 6_identify_ptu.py script"
    # Call 6_identify_ptu.py script
    python3 6_identify_ptu.py "$gff" possible_ptu.txt
    echo "Identify PTU regions"

    echo "Running 7_parser_UTR_sense_antisense.py script"
    # Call 7_parser_UTR_sense_antisense.py script
    python3 7_parser_UTR_sense_antisense.py "${output_folder}/annotation_ncRNA_final.txt" "${output_folder}/annotation_lncRNA.txt" possible_ptu.txt "${output_folder}/ncRNAs_location_direction.bed" "$utr5" "$utr3"
    echo "Annotation transcript at sense and location level"

    echo "Running bedtools getfasta"
    #Extract fasta sequences 
    bedtools getfasta -fi "$fasta" -bed "${output_folder}/ncRNAs_location_direction.bed" -fo "${output_folder}/all_ncrna.fasta" -name+
    echo "Extracted fasta sequences"

    #Running filter blastx
    echo "Running diamond"
    diamond blastx --query "${output_folder}/all_ncrna.fasta" --db "${database}/pfam_database.dmnd" --out "${output_folder}/ncrna_pfam-cov80-max1.tab" --outfmt 6 qseqid sseqid pident qcovhsp length qlen slen qstart qend sstart send evalue bitscore stitle --id 90 --query-cover 80 --evalue 1e-5 --threads "$threads" --max-target-seqs 1
    echo "diamond blastx done successefully"
    
    echo "Final selection of ncRNA"
    python3 8_select_ncrna.py "${output_folder}/ncRNAs_location_direction.bed" "${output_folder}/ncrna_pfam-cov80-max1.tab" "${output_folder}/final_ncRNAs_location_direction.bed"

done < "$dir_list"

echo "cat all bed outputs, sort and merge all unique anotates ncRNA"
allstages=""
# Loop through directories ending with "_out"
for subdir in *_out; do
    # Ensure it's a directory
    if [ -d "$subdir" ]; then
        # Check if the expected file exists
        if [ -f "$subdir/final_ncRNAs_location_direction.bed" ]; then
            allstages+="$subdir/final_ncRNAs_location_direction.bed "
        else
            echo "Warning: File not found in $subdir: final_ncRNAs_location_direction.bed"
        fi
    fi
done
echo "${allstages}" #Print the concatenated result


cat ${allstages} > "${output_base}/final_all_ncrna.bed"
sort -k1,1 -k2,2n "${output_base}/final_all_ncrna.bed" > "${output_base}/sorted_all_ncRNA.bed"
bedtools merge -i "${output_base}/sorted_all_ncRNA.bed" -s -c 4,6,7,8,9,10 -o distinct,distinct,distinct,distinct,distinct,distinct > "${output_base}/unique_sort_allncrna.tab"

python3 9_remake_output.py "${output_base}/unique_sort_allncrna.tab" "${output_base}/df_allncrna.tab" "${output_base}/df_allncrna.bed" "${output_base}/df_allncrna.gff"

echo "Extract fasta non-coding RNA"
#Extract fasta sequences 
bedtools getfasta -fi "$fasta" -bed "${output_base}/df_allncrna.bed" -fo "${output_base}/fasta_ncrna.fasta" -name+
echo "Extracted fasta sequences"