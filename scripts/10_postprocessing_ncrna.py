import sys
import pandas as pd
import re

# python3 10_postprocessing_ncrna.py df_allncrna.tab fasta_ncrna.fasta_results_all.scores output_ptrnapred1.txt tRNAscan-output.tab output_snoscan.txt output_rnacon.txt df_allncrna_final.tab
def process_file(ncrna_tab, out_portrait, out_tpredrna, out_trnascan, out_snoscan, out_rnacon, out_ncrna_tab):

    # Importar o arquivo usando pandas
    ncrna = pd.read_csv(ncrna_tab, header=None,  sep='\t', names=["nc_chr","nc_coordi","nc_coordf","match_id","nc_strand","match_length","nc_sentido","nc_position", "nc_lenght", "nc_id"])
    #post-processing the output PORTRAIT
    output = pd.read_csv(out_portrait, names = ["column"])
    #  extração deve começar imediatamente ">" ate ":" segundo grupo captura o valor decimal após os dois pontos,
    output[['nc_id', 'score']] = output['column'].str.extract(r'>([\w_]+).*:.*? (\d\.\d+)$')
    #output['score'] = output['score'].astype(float)
    merged_df = pd.merge(ncrna, output[['nc_id', 'score']], on='nc_id', how='left')
    
    #post-processing the output tRNApred
    output_tpred = pd.read_csv(out_tpredrna, names = ["column"])
    selec1 = output_tpred[output_tpred['column'].str.contains(r'^(SEQ NAME:>)', regex=True)]
    selec2 = output_tpred[output_tpred['column'].str.contains(r'^(RNA Family:)', regex=True)]
    rnaid = selec1['column'].str.extract(r'SEQ NAME:>([\w_]+)(?=::)')  
    rnafamily = selec2['column'].str.extract(r'RNA Family:([\w_]+)')
    rnaid.reset_index(drop=True, inplace=True)
    rnafamily.reset_index(drop=True, inplace=True)
    output_tpred = pd.concat([rnaid, rnafamily], axis=1)
    output_tpred.columns = ['nc_id', 'rna_family']
    merged_df2 = pd.merge(merged_df, output_tpred, on='nc_id', how='left')
    #post-processing the output tRNAscan
    trna = pd.read_csv(out_trnascan, sep='\t')
    trna = trna.iloc[2:,:1]
    trna = trna.reset_index()
    trna.columns= ["name", "tRNAscan"]
    trna[['nc_id']] = trna['name'].str.extract(r'^([\w_.]+)::')
    trna_uniq = trna.drop_duplicates(subset=['nc_id'])                
    merged_df3 = pd.merge(merged_df2, trna_uniq[['nc_id', 'tRNAscan']], on='nc_id', how='left')
    merged_df3['tRNAscan'] = merged_df3['tRNAscan'].fillna(0)
    
    #post-processing the output snoscan 
    output_snoscan = pd.read_csv(out_snoscan, names=["column"])
    snoscan = output_snoscan[output_snoscan['column'].str.contains(">>")].copy()
    snoscan[['nc_id', 'snoscan']] = snoscan['column'].str.extract(r'>>\s*([\w\d_]+)::.*?\s+([\d]+\.\d+)')
    snoscan.reset_index(drop=True, inplace=True)
    snoscan_uniq = snoscan.drop_duplicates(subset=['nc_id'])
    merged_df4 = pd.merge(merged_df3, snoscan_uniq[['nc_id', 'snoscan']], on='nc_id', how='left')
    merged_df4['snoscan'] = merged_df4['snoscan'].fillna(0)

    #post=processing the output rnacon

    with open(out_rnacon, encoding='utf-8') as f:
        lines = f.readlines()[7:]
    lines = [line.strip() for line in lines if line.strip() and not line.startswith('#')]

    # Processa cada linha separando por dois ou mais espaços ou tabulações
    data = []
    for line in lines:
        parts = re.split(r'\s{2,}|\t+', line)
        while len(parts) < 3:
            parts.append(None)
        data.append(parts[:3])  # Garante no máximo 3 colunas

    # Cria DataFrame
    rnacon = pd.DataFrame(data, columns=["id", "rnacon", "Classification_rnacon"])
    rnacon['nc_id'] = rnacon['id'].str.extract(r'\.(.*?)::')
    rnacon['nc_id'] = rnacon['nc_id'].str.strip().astype(str)
    merged_df5 = pd.merge(merged_df4, rnacon[['nc_id', 'rnacon']], on='nc_id', how='left')
    merged_df5.to_csv(out_ncrna_tab, sep='\t', index=False, header=False)

# Importar e executar o script
if __name__ == "__main__":
    ncrna_tab = sys.argv[1]
    out_portrait = sys.argv[2]
    out_tpredrna = sys.argv[3]
    out_trnascan = sys.argv[4]
    out_snoscan = sys.argv[5]
    out_rnacon = sys.argv[6]
    out_ncrna_tab = sys.argv[7]
    process_file(ncrna_tab, out_portrait, out_tpredrna, out_trnascan, out_snoscan, out_rnacon, out_ncrna_tab)

