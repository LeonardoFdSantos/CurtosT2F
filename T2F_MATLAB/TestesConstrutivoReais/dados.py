import numpy as np
import pandas as pd

dadosFimLinha = pd.read_csv('ResultadoFimLinha.csv')
dadosMeioLinha1 = pd.read_csv('ResultadoMeioLinha1.csv')
dadosMeioLinha2 = pd.read_csv('ResultadoMeioLinha2.csv')

dadosFimLinha_SemFuse = pd.read_csv('ResultadoFimLinha_SemFuse.csv')
dadosMeioLinha1_SemFuse = pd.read_csv('ResultadoMeioLinha1_SemFuse.csv')
dadosMeioLinha2_SemFuse = pd.read_csv('ResultadoMeioLinha2_SemFuse.csv')

dadosFimLinha = dadosFimLinha.replace(np.inf, 0)
dadosMeioLinha1 = dadosMeioLinha1.replace(np.inf, 0)
dadosMeioLinha2 = dadosMeioLinha2.replace(np.inf, 0)

dadosFimLinha_SemFuse = dadosFimLinha_SemFuse.replace(np.inf, 0)
dadosMeioLinha1_SemFuse = dadosMeioLinha1_SemFuse.replace(np.inf, 0)
dadosMeioLinha2_SemFuse = dadosMeioLinha2_SemFuse.replace(np.inf, 0)

with pd.ExcelWriter('Dados.xlsx') as writer:  
    dadosFimLinha.to_excel(writer, sheet_name='Fim Linha', index=False)
    dadosMeioLinha1.to_excel(writer, sheet_name='Meio Linha 1', index=False)
    dadosMeioLinha2.to_excel(writer, sheet_name='Meio Linha 2', index=False)

with pd.ExcelWriter('Dados_SemFuse.xlsx') as writer:  
    dadosFimLinha_SemFuse.to_excel(writer, sheet_name='Fim Linha Sem Fuse', index=False)
    dadosMeioLinha1_SemFuse.to_excel(writer, sheet_name='Meio Linha 1 Sem Fuse', index=False)
    dadosMeioLinha2_SemFuse.to_excel(writer, sheet_name='Meio Linha 2 Sem Fuse', index=False)
'''
dadosGhendyanoT2F = pd.read_csv('TesteEloGhendy_16.csv')
dadosGhendyanoT2FConsumidor = pd.read_csv('TesteEloGhendyConsumidor_16.csv')

dadosGhendyanoT2F = dadosGhendyanoT2F.replace(np.inf, 0)
dadosGhendyanoT2FConsumidor = dadosGhendyanoT2FConsumidor.replace(np.inf, 0)

with pd.ExcelWriter('DadosFuse.xlsx') as writer:
    dadosGhendyanoT2F.to_excel(writer, sheet_name='Linha', index=False)
    dadosGhendyanoT2FConsumidor.to_excel(writer, sheet_name='Consumidor', index=False)
'''