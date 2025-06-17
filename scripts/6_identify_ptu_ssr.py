import sys

# python 6_identify_ptu.py TriTrypDB-30_LbraziliensisMHOMBR75M2903.gff possible_ptu.txt possible_ssr.txt
chr_old = ""
start_old = ""
end_old = ""
strand_old = ""

def process_file(gff, out_file, out_file2): #possible_ptu.txt possible_ssr.txt are output files
    global chr_old, start_old, end_old, strand_old

    with open(gff, 'r') as fh:
        # Process header gff and save chromossomes and start and end position
        chr_bounds = {}
        cds = []
        for line in fh:
            if line.startswith("##sequence-region"):
                fields = line.strip().split()
                chr_bounds[fields[1]] = (int(fields[2]), int(fields[3]))  # {chr: (start_gff, end_gff)}
            elif not line.startswith("##"):
                cds.append(line.strip())  # Save rows that not to be comments

    # Sort rows by chromossomal (column 0) and start positions (column third) 
    cds.sort(key=lambda row: (row.split("\t")[0], int(row.split("\t")[3])))

    with open(out_file, 'w') as ptu, open(out_file2, 'w') as ssr:
        
        for row in cds:
            fields = row.split("\t")
            if len(fields) < 7:
                continue
            chr, mol, start, end, strand = fields[0], fields[2], int(fields[3]), int(fields[4]), fields[6]
            if mol != "CDS":  # Continue with CDS rows
                continue
            # Update start and end position from old or new chromossome
            if chr != chr_old:  #New chromossome
                if chr_old in chr_bounds:
                    final_chr_pos = chr_bounds[chr_old][1]
                    ssr.write(f"{chr_old}\t.\t{end_old + 1}\t{final_chr_pos}\tChromossomal end\n") #Save tha last region of chromossome   
                if chr_old and strand_old:  # Save information for the previous chromosome
                    ptu.write(f"{chr_old}\t{strand_old}\t{start_old}\t{end_old}\tPTU\n") ###COMMENT THIS LINE
                if start > 1: #Save the position start in the change of chromossome
                    ssr.write(f"{chr}\t.\t1\t{start - 1}\tChromossomal start\n") 
                # Atualizar valores para o novo cromossomo
                chr_old = chr
                start_old = start  # Sempre 1 na primeira vez
                end_old = end
                strand_old = strand
            elif strand != strand_old:  # Change of strand in the same chromossome
                ptu.write(f"{chr_old}\t{strand_old}\t{start_old}\t{end_old}\tPTU\n")
                end_old = end_old + 1 #position of last CDs
                start = start -1 #position of next CDs
                if strand_old == "+": #Give information for convergence of PTUs
                    ssr.write(f"{chr_old}\t.\t{end_old}\t{start}\tcSSR\n")
                if strand_old == "-":
                    ssr.write(f"{chr_old}\t.\t{end_old}\t{start}\tdSSR\n")
                start_old = start  #Update start
                strand_old = strand
                end_old = end
            else:
                end_old = end  # Update final position
                
        # Save the last line if valid
        if chr_old and strand_old:
            ptu.write(f"{chr_old}\t{strand_old}\t{start_old}\t{end_old}\tPTU\n") 

# Importar e executar o script
if __name__ == "__main__":
    gff = sys.argv[1]
    out_file = sys.argv[2]
    out_file2 = sys.argv[3]
    process_file(gff, out_file, out_file2)
