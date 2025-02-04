import sys
import pandas as pd

# python3 10_postprocessing_ncrna.py df_allncrna.tab fasta_ncrna.fasta_results_all.scores output_ptrnapred1.txt tRNAscan-output.tab
def process_file(ncrna_tab, out_portrait, out_tpredrna, out_trnascan):

    # Importar o arquivo usando pandas
    ncrna = pd.read_csv(ncrna_tab, header=None,  sep='\t', names=["nc_chr","nc_coordi","nc_coordf","match_id","nc_strand","match_length","nc_location","nc_sentido","nc_position", "nc_lenght", "nc_id"])
    output = pd.read_csv(out_portrait, names = ["column"])
    #  extração deve começar imediatamente ">" ate ":" segundo grupo captura o valor decimal após os dois pontos,
    output[['nc_id', 'score']] = output['column'].str.extract(r'>([\w_]+).*:.*? (\d\.\d+)$')
    #output['score'] = output['score'].astype(float)

    merged_df = pd.merge(ncrna, output[['nc_id', 'score']], on='nc_id', how='left')
    

    output_tpred = pd.read_csv(out_tpredrna, names = ["column"])
    selec1 = output_tpred[output_tpred['column'].str.contains(r'^(SEQ NAME:>)', regex=True)]
    selec2 = output_tpred[output_tpred['column'].str.contains(r'^(RNA Family:)', regex=True)]
    rnaid = selec1['column'].str.extract(r'SEQ NAME:>([\w_]+)(?=::)')  
    rnafamily = selec2['column'].str.extract(r'RNA Family:([\w_]+)')
    rnaid.reset_index(drop=True, inplace=True)
    rnafamily.reset_index(drop=True, inplace=True)
    # Concatenar rnaid e rnafamily de forma alinhada
    output_tpred = pd.concat([rnaid, rnafamily], axis=1)

    # Definir os nomes das colunas para o DataFrame
    output_tpred.columns = ['nc_id', 'rna_family']
    merged_df2 = pd.merge(merged_df, output_tpred, on='nc_id', how='left')

    trna = pd.read_csv(out_trnascan, sep='\t')
    trna = trna.iloc[2:,:1]
    trna = trna.reset_index()
    # Definir a segunda linha como cabeçalho
    trna.columns= ["name", "tRNAscan"]
    trna[['nc_id']] = trna['name'].str.extract(r'^([\w_.]+)::')
                            
    merged_df3 = pd.merge(merged_df2, trna[['nc_id', 'tRNAscan']], on='nc_id', how='left')
    merged_df3['tRNAscan'] = merged_df3['tRNAscan'].fillna(0)
    merged_df3.to_csv(ncrna_tab, sep='\t', index=False, header=False)
    print(merged_df3)
# Importar e executar o script
if __name__ == "__main__":
    ncrna_tab = sys.argv[1]
    out_portrait = sys.argv[2]
    out_tpredrna = sys.argv[3]
    out_trnascan = sys.argv[4]
    process_file(ncrna_tab, out_portrait, out_tpredrna, out_trnascan)