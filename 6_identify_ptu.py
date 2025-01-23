import sys

# python 5_identify_ptu.py TriTrypDB-30_LbraziliensisMHOMBR75M2903.gff possible_ptu.txt
chr_old = ""
start_old = ""
end_old = ""
strand_old = ""

def process_file(gff, out_file):
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
                cds.append(line.strip())  # Save lines that not to be comments

    # Sort rows by chromossomal (column 0) and start positions (column third) 
    cds.sort(key=lambda row: (row.split("\t")[0], int(row.split("\t")[3])))

    with open(out_file, 'w') as ptu:
        for row in cds:
            fields = row.split("\t")
            
            if len(fields) < 7:
                continue
            chr, mol, start, end, strand = fields[0], fields[2], int(fields[3]), int(fields[4]), fields[6]
            
            # Conitue with CDS lines
            if mol != "CDS":
                continue
            
            # Atualizar a posição inicial e final com base no cromossomo
            if chr != chr_old:  # Novo cromossomo
                if chr_old and strand_old:  # Save information for the previous chromosome
                    end_old = chr_bounds.get(chr_old, (None, end_old))[1]  # Use the chromosome bounds if available
                    ptu.write(f"{chr_old}\t{strand_old}\t{start_old}\t{end_old}\n")
                # Atualizar valores para o novo cromossomo
                chr_old = chr
                start_old = 1  # Sempre 1 na primeira vez
                end_old = end
                strand_old = strand
            elif strand != strand_old:  # Mudança de strand no mesmo cromossomo
                end_old = start - 1  # Atualizar a posição final para a mudança de strand
                ptu.write(f"{chr_old}\t{strand_old}\t{start_old}\t{end_old}\n")
                start_old = start  # Atualizar o início para o próximo segmento
                strand_old = strand
            else:
                end_old = end  # Atualizar a posição final no mesmo segmento
        
        # Save the last line if valid
        if chr_old and strand_old:
            end_old = chr_bounds.get(chr_old, (None, end_old))[1]  # Use chromosome bounds if available
            ptu.write(f"{chr_old}\t{strand_old}\t{start_old}\t{end_old}\n")


# Importar e executar o script
if __name__ == "__main__":
    gff = sys.argv[1]
    out_file = sys.argv[2]
    process_file(gff, out_file)
