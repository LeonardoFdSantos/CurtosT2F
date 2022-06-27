import pandas as pd
import numpy as np

valores_resultados_calculo_fim_linha_com_compensacao = pd.read_csv('valores_resultados_calculo_fim_linha_com_compensacao.csv');
valores_resultados_calculo_fim_linha_sem_compensacao = pd.read_csv('valores_resultados_calculo_fim_linha_sem_compensacao.csv');
valores_resultados_calculo_meio_linha_com_compensacao = pd.read_csv('valores_resultados_calculo_meio_linha_com_compensacao.csv');
valores_resultados_calculo_meio_linha_sem_compensacao = pd.read_csv('valores_resultados_calculo_meio_linha_sem_compensacao.csv');

valores_resultados_simulacao_fim_linha_com_compensacao = pd.read_csv('valores_resultados_simulacao_fim_linha_com_compensacao.csv');
valores_resultados_simulacao_fim_linha_sem_compensacao = pd.read_csv('valores_resultados_simulacao_fim_linha_sem_compensacao.csv');
valores_resultados_simulacao_meio_linha_com_compensacao = pd.read_csv('valores_resultados_simulacao_meio_linha_com_compensacao.csv');
valores_resultados_simulacao_meio_linha_sem_compensacao = pd.read_csv('valores_resultados_simulacao_meio_linha_sem_compensacao.csv');

valores_resultados_calculo_fim_linha_com_compensacao[''] = '';
valores_resultados_calculo_fim_linha_sem_compensacao[''] = '';
valores_resultados_calculo_meio_linha_com_compensacao[''] = '';
valores_resultados_calculo_meio_linha_sem_compensacao[''] = '';

df_resultados_fim_linha_com_compensacao = pd.concat([valores_resultados_calculo_fim_linha_com_compensacao, valores_resultados_simulacao_fim_linha_com_compensacao], axis=1);
df_resultados_fim_linha_sem_compensacao = pd.concat([valores_resultados_calculo_fim_linha_sem_compensacao, valores_resultados_simulacao_fim_linha_sem_compensacao], axis=1);
df_resultados_meio_linha_com_compensacao = pd.concat([valores_resultados_calculo_meio_linha_com_compensacao, valores_resultados_simulacao_meio_linha_com_compensacao], axis=1);
df_resultados_meio_linha_sem_compensacao = pd.concat([valores_resultados_calculo_meio_linha_sem_compensacao, valores_resultados_simulacao_meio_linha_sem_compensacao], axis=1);

localizacao = 'DadosT2F.xlsx'

with pd.ExcelWriter(localizacao) as writer:
    df_resultados_fim_linha_com_compensacao.to_excel(writer, sheet_name='Fim_Linha_com_Comp');
    df_resultados_fim_linha_sem_compensacao.to_excel(writer, sheet_name='Fim_Linha_sem_Comp');
    df_resultados_meio_linha_com_compensacao.to_excel(writer, sheet_name='Meio_Linha_com_Comp');
    df_resultados_meio_linha_sem_compensacao.to_excel(writer, sheet_name='Meio_Linha_sem_Comp');

print('Dados exportados, para '+ localizacao);
print('Terminado!');