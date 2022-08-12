import pandas as pd
import numpy as np

valores_resultados_fim_linha_com_compensacao = pd.read_csv('valores_resultados_fim_linha_com_compensacao.csv');
valores_resultados_meio_linha_com_compensacao = pd.read_csv('valores_resultados_meio_linha_com_compensacao.csv');
valores_resultados_meio_linha_sem_compensacao = pd.read_csv('valores_resultados_meio_linha_sem_compensacao.csv');
valores_resultados_fim_linha_sem_compensacao = pd.read_csv('valores_resultados_fim_linha_sem_compensacao.csv');

valores_resultados_simulacao_fim_linha = pd.read_csv('valores_resultados_simulacao_fim_linha_sem_comp.csv');
valores_resultados_simulacao_meio_linha = pd.read_csv('valores_resultados_simulacao_meio_linha_sem_comp.csv');
valores_resultados_simulacao_meio_linha_com_comp = pd.read_csv('valores_resultados_simulacao_meio_linha_com_comp.csv');
valores_resultados_simulacao_fim_linha_com_comp = pd.read_csv('valores_resultados_simulacao_fim_linha_com_comp.csv');

valores_resultados_meio_linha_sem_compensacao['Erro_IA'] = (valores_resultados_meio_linha_sem_compensacao['IA_Sim'] - valores_resultados_meio_linha_sem_compensacao['IA_Calc'])/valores_resultados_meio_linha_sem_compensacao['IA_Calc']
valores_resultados_meio_linha_sem_compensacao['Erro_IB'] = (valores_resultados_meio_linha_sem_compensacao['IB_Sim'] - valores_resultados_meio_linha_sem_compensacao['IB_Calc'])/valores_resultados_meio_linha_sem_compensacao['IB_Calc']
valores_resultados_meio_linha_sem_compensacao['Erro_IC'] = (valores_resultados_meio_linha_sem_compensacao['IC_Sim'] - valores_resultados_meio_linha_sem_compensacao['IC_Calc'])/valores_resultados_meio_linha_sem_compensacao['IC_Calc']

valores_resultados_fim_linha_sem_compensacao['Erro_IA'] = (valores_resultados_fim_linha_sem_compensacao['IA_Sim'] - valores_resultados_fim_linha_sem_compensacao['IA_Calc'])/valores_resultados_fim_linha_sem_compensacao['IA_Calc']
valores_resultados_fim_linha_sem_compensacao['Erro_IB'] = (valores_resultados_fim_linha_sem_compensacao['IB_Sim'] - valores_resultados_fim_linha_sem_compensacao['IB_Calc'])/valores_resultados_fim_linha_sem_compensacao['IB_Calc']
valores_resultados_fim_linha_sem_compensacao['Erro_IC'] = (valores_resultados_fim_linha_sem_compensacao['IC_Sim'] - valores_resultados_fim_linha_sem_compensacao['IC_Calc'])/valores_resultados_fim_linha_sem_compensacao['IC_Calc']

valores_resultados_fim_linha_com_compensacao['Erro_IA'] = (valores_resultados_fim_linha_com_compensacao['IA_Sim'] - valores_resultados_fim_linha_com_compensacao['IA_Calc'])/valores_resultados_fim_linha_com_compensacao['IA_Calc']
valores_resultados_fim_linha_com_compensacao['Erro_IB'] = (valores_resultados_fim_linha_com_compensacao['IB_Sim'] - valores_resultados_fim_linha_com_compensacao['IB_Calc'])/valores_resultados_fim_linha_com_compensacao['IB_Calc']
valores_resultados_fim_linha_com_compensacao['Erro_IC'] = (valores_resultados_fim_linha_com_compensacao['IC_Sim'] - valores_resultados_fim_linha_com_compensacao['IC_Calc'])/valores_resultados_fim_linha_com_compensacao['IC_Calc']

valores_resultados_meio_linha_com_compensacao['Erro_IA'] = (valores_resultados_meio_linha_com_compensacao['IA_Sim'] - valores_resultados_meio_linha_com_compensacao['IA_Calc'])/valores_resultados_meio_linha_com_compensacao['IA_Calc']
valores_resultados_meio_linha_com_compensacao['Erro_IB'] = (valores_resultados_meio_linha_com_compensacao['IB_Sim'] - valores_resultados_meio_linha_com_compensacao['IB_Calc'])/valores_resultados_meio_linha_com_compensacao['IB_Calc']
valores_resultados_meio_linha_com_compensacao['Erro_IC'] = (valores_resultados_meio_linha_com_compensacao['IC_Sim'] - valores_resultados_meio_linha_com_compensacao['IC_Calc'])/valores_resultados_meio_linha_com_compensacao['IC_Calc']

valores_resultados_simulacao_fim_linha['Erro_IA'] = (valores_resultados_simulacao_fim_linha['IA_T2F'] - valores_resultados_simulacao_fim_linha['IA_TRIF'])/valores_resultados_simulacao_fim_linha['IA_TRIF']
valores_resultados_simulacao_fim_linha['Erro_IB'] = (valores_resultados_simulacao_fim_linha['IB_T2F'] - valores_resultados_simulacao_fim_linha['IB_TRIF'])/valores_resultados_simulacao_fim_linha['IB_TRIF']
valores_resultados_simulacao_fim_linha['Erro_IC'] = (valores_resultados_simulacao_fim_linha['IC_T2F'] - valores_resultados_simulacao_fim_linha['IC_TRIF'])/valores_resultados_simulacao_fim_linha['IC_TRIF']

valores_resultados_simulacao_fim_linha_com_comp['Erro_IA'] = (valores_resultados_simulacao_fim_linha_com_comp['IA_T2F'] - valores_resultados_simulacao_fim_linha_com_comp['IA_TRIF'])/valores_resultados_simulacao_fim_linha_com_comp['IA_TRIF']
valores_resultados_simulacao_fim_linha_com_comp['Erro_IB'] = (valores_resultados_simulacao_fim_linha_com_comp['IB_T2F'] - valores_resultados_simulacao_fim_linha_com_comp['IB_TRIF'])/valores_resultados_simulacao_fim_linha_com_comp['IB_TRIF']
valores_resultados_simulacao_fim_linha_com_comp['Erro_IC'] = (valores_resultados_simulacao_fim_linha_com_comp['IC_T2F'] - valores_resultados_simulacao_fim_linha_com_comp['IC_TRIF'])/valores_resultados_simulacao_fim_linha_com_comp['IC_TRIF']

valores_resultados_simulacao_meio_linha['Erro_IA'] = (valores_resultados_simulacao_meio_linha['IA_T2F'] - valores_resultados_simulacao_meio_linha['IA_TRIF'])/valores_resultados_simulacao_meio_linha['IA_TRIF']
valores_resultados_simulacao_meio_linha['Erro_IB'] = (valores_resultados_simulacao_meio_linha['IB_T2F'] - valores_resultados_simulacao_meio_linha['IB_TRIF'])/valores_resultados_simulacao_meio_linha['IB_TRIF']
valores_resultados_simulacao_meio_linha['Erro_IC'] = (valores_resultados_simulacao_meio_linha['IC_T2F'] - valores_resultados_simulacao_meio_linha['IC_TRIF'])/valores_resultados_simulacao_meio_linha['IC_TRIF']

valores_resultados_simulacao_meio_linha_com_comp['Erro_IA'] = (valores_resultados_simulacao_meio_linha_com_comp['IA_T2F'] - valores_resultados_simulacao_meio_linha_com_comp['IA_TRIF'])/valores_resultados_simulacao_meio_linha_com_comp['IA_TRIF']
valores_resultados_simulacao_meio_linha_com_comp['Erro_IB'] = (valores_resultados_simulacao_meio_linha_com_comp['IB_T2F'] - valores_resultados_simulacao_meio_linha_com_comp['IB_TRIF'])/valores_resultados_simulacao_meio_linha_com_comp['IB_TRIF']
valores_resultados_simulacao_meio_linha_com_comp['Erro_IC'] = (valores_resultados_simulacao_meio_linha_com_comp['IC_T2F'] - valores_resultados_simulacao_meio_linha_com_comp['IC_TRIF'])/valores_resultados_simulacao_meio_linha_com_comp['IC_TRIF']

valores_resultados_meio_linha_sem_compensacao_Trif = valores_resultados_meio_linha_sem_compensacao.loc[valores_resultados_meio_linha_sem_compensacao['tipoCurto'] == 1]
valores_resultados_meio_linha_sem_compensacao_AC = valores_resultados_meio_linha_sem_compensacao.loc[valores_resultados_meio_linha_sem_compensacao['tipoCurto'] == 2]
valores_resultados_meio_linha_sem_compensacao_BC = valores_resultados_meio_linha_sem_compensacao.loc[valores_resultados_meio_linha_sem_compensacao['tipoCurto'] == 3]
valores_resultados_meio_linha_sem_compensacao_AB = valores_resultados_meio_linha_sem_compensacao.loc[valores_resultados_meio_linha_sem_compensacao['tipoCurto'] == 4]

valores_resultados_meio_linha_com_compensacao_Trif = valores_resultados_meio_linha_com_compensacao.loc[valores_resultados_meio_linha_com_compensacao['tipoCurto'] == 1]
valores_resultados_meio_linha_com_compensacao_AC = valores_resultados_meio_linha_com_compensacao.loc[valores_resultados_meio_linha_com_compensacao['tipoCurto'] == 2]
valores_resultados_meio_linha_com_compensacao_BC = valores_resultados_meio_linha_com_compensacao.loc[valores_resultados_meio_linha_com_compensacao['tipoCurto'] == 3]
valores_resultados_meio_linha_com_compensacao_AB = valores_resultados_meio_linha_com_compensacao.loc[valores_resultados_meio_linha_com_compensacao['tipoCurto'] == 4]

valores_resultados_simulacao_meio_linha_Trif = valores_resultados_simulacao_meio_linha.loc[valores_resultados_simulacao_meio_linha['tipoCurto'] == 1]
valores_resultados_simulacao_meio_linha_AC = valores_resultados_simulacao_meio_linha.loc[valores_resultados_simulacao_meio_linha['tipoCurto'] == 2]
valores_resultados_simulacao_meio_linha_BC = valores_resultados_simulacao_meio_linha.loc[valores_resultados_simulacao_meio_linha['tipoCurto'] == 3]
valores_resultados_simulacao_meio_linha_AB = valores_resultados_simulacao_meio_linha.loc[valores_resultados_simulacao_meio_linha['tipoCurto'] == 4]
valores_resultados_meio_linha_sem_curto = valores_resultados_simulacao_meio_linha.loc[valores_resultados_simulacao_meio_linha['tipoCurto'] == 5]

valores_resultados_simulacao_meio_linha_Trif_com_comp = valores_resultados_simulacao_meio_linha_com_comp.loc[valores_resultados_simulacao_meio_linha_com_comp['tipoCurto'] == 1]
valores_resultados_simulacao_meio_linha_AC_com_comp = valores_resultados_simulacao_meio_linha_com_comp.loc[valores_resultados_simulacao_meio_linha_com_comp['tipoCurto'] == 2]
valores_resultados_simulacao_meio_linha_BC_com_comp = valores_resultados_simulacao_meio_linha_com_comp.loc[valores_resultados_simulacao_meio_linha_com_comp['tipoCurto'] == 3]
valores_resultados_simulacao_meio_linha_AB_com_comp = valores_resultados_simulacao_meio_linha_com_comp.loc[valores_resultados_simulacao_meio_linha_com_comp['tipoCurto'] == 4]
valores_resultados_meio_linha_sem_curto_com_comp = valores_resultados_simulacao_meio_linha_com_comp.loc[valores_resultados_simulacao_meio_linha_com_comp['tipoCurto'] == 5]

valores_resultados_meio_linha_sem_compensacao_Trif = valores_resultados_meio_linha_sem_compensacao_Trif.drop(['tipoCurto'], axis=1)
valores_resultados_meio_linha_sem_compensacao_AC = valores_resultados_meio_linha_sem_compensacao_AC.drop(['tipoCurto'], axis=1)
valores_resultados_meio_linha_sem_compensacao_BC = valores_resultados_meio_linha_sem_compensacao_BC.drop(['tipoCurto'], axis=1)
valores_resultados_meio_linha_sem_compensacao_AB = valores_resultados_meio_linha_sem_compensacao_AB.drop(['tipoCurto'], axis=1)

valores_resultados_meio_linha_com_compensacao_Trif = valores_resultados_meio_linha_com_compensacao_Trif.drop(['tipoCurto'], axis=1)
valores_resultados_meio_linha_com_compensacao_AC = valores_resultados_meio_linha_com_compensacao_AC.drop(['tipoCurto'], axis=1)
valores_resultados_meio_linha_com_compensacao_BC = valores_resultados_meio_linha_com_compensacao_BC.drop(['tipoCurto'], axis=1)
valores_resultados_meio_linha_com_compensacao_AB = valores_resultados_meio_linha_com_compensacao_AB.drop(['tipoCurto'], axis=1)

valores_resultados_simulacao_meio_linha_Trif = valores_resultados_simulacao_meio_linha_Trif.drop(['tipoCurto'], axis=1)
valores_resultados_simulacao_meio_linha_AC = valores_resultados_simulacao_meio_linha_AC.drop(['tipoCurto'], axis=1)
valores_resultados_simulacao_meio_linha_BC = valores_resultados_simulacao_meio_linha_BC.drop(['tipoCurto'], axis=1)
valores_resultados_simulacao_meio_linha_AB = valores_resultados_simulacao_meio_linha_AB.drop(['tipoCurto'], axis=1)
valores_resultados_meio_linha_sem_curto = valores_resultados_meio_linha_sem_curto.drop(['tipoCurto'], axis=1)

valores_resultados_simulacao_meio_linha_Trif_com_comp = valores_resultados_simulacao_meio_linha_Trif_com_comp.drop(['tipoCurto'], axis=1)
valores_resultados_simulacao_meio_linha_AC_com_comp = valores_resultados_simulacao_meio_linha_AC_com_comp.drop(['tipoCurto'], axis=1)
valores_resultados_simulacao_meio_linha_BC_com_comp = valores_resultados_simulacao_meio_linha_BC_com_comp.drop(['tipoCurto'], axis=1)
valores_resultados_simulacao_meio_linha_AB_com_comp = valores_resultados_simulacao_meio_linha_AB_com_comp.drop(['tipoCurto'], axis=1)
valores_resultados_meio_linha_sem_curto_com_comp = valores_resultados_meio_linha_sem_curto_com_comp.drop(['tipoCurto'], axis=1)

valores_resultados_fim_linha_sem_compensacao = valores_resultados_fim_linha_sem_compensacao.set_index(['tipoCurto'])
valores_resultados_fim_linha_com_compensacao = valores_resultados_fim_linha_com_compensacao.set_index(['tipoCurto'])
valores_resultados_simulacao_fim_linha = valores_resultados_simulacao_fim_linha.set_index(['tipoCurto'])

valores_resultados_meio_linha_sem_compensacao_Trif = valores_resultados_meio_linha_sem_compensacao_Trif.set_index(['n'])
valores_resultados_meio_linha_sem_compensacao_AC = valores_resultados_meio_linha_sem_compensacao_AC.set_index(['n'])
valores_resultados_meio_linha_sem_compensacao_BC = valores_resultados_meio_linha_sem_compensacao_BC.set_index(['n'])
valores_resultados_meio_linha_sem_compensacao_AB = valores_resultados_meio_linha_sem_compensacao_AB.set_index(['n'])

valores_resultados_meio_linha_com_compensacao_Trif = valores_resultados_meio_linha_com_compensacao_Trif.set_index(['n'])
valores_resultados_meio_linha_com_compensacao_AC = valores_resultados_meio_linha_com_compensacao_AC.set_index(['n'])
valores_resultados_meio_linha_com_compensacao_BC = valores_resultados_meio_linha_com_compensacao_BC.set_index(['n'])
valores_resultados_meio_linha_com_compensacao_AB = valores_resultados_meio_linha_com_compensacao_AB.set_index(['n'])

valores_resultados_simulacao_meio_linha_Trif = valores_resultados_simulacao_meio_linha_Trif.set_index(['n'])
valores_resultados_simulacao_meio_linha_AC = valores_resultados_simulacao_meio_linha_AC.set_index(['n'])
valores_resultados_simulacao_meio_linha_BC = valores_resultados_simulacao_meio_linha_BC.set_index(['n'])
valores_resultados_simulacao_meio_linha_AB = valores_resultados_simulacao_meio_linha_AB.set_index(['n'])
valores_resultados_meio_linha_sem_curto = valores_resultados_meio_linha_sem_curto.set_index(['n'])

valores_resultados_simulacao_meio_linha_Trif_com_comp = valores_resultados_simulacao_meio_linha_Trif_com_comp.set_index(['n'])
valores_resultados_simulacao_meio_linha_AC_com_comp = valores_resultados_simulacao_meio_linha_AC_com_comp.set_index(['n'])
valores_resultados_simulacao_meio_linha_BC_com_comp = valores_resultados_simulacao_meio_linha_BC_com_comp.set_index(['n'])
valores_resultados_simulacao_meio_linha_AB_com_comp = valores_resultados_simulacao_meio_linha_AB_com_comp.set_index(['n'])
valores_resultados_meio_linha_sem_curto_com_comp = valores_resultados_meio_linha_sem_curto_com_comp.set_index(['n'])

localizacao_T2F = 'DadosT2F.xlsx'
with pd.ExcelWriter(localizacao_T2F) as writer:
    valores_resultados_fim_linha_sem_compensacao.to_excel(writer, sheet_name='Fim_Linha_sem_Comp');    
    valores_resultados_meio_linha_sem_compensacao_Trif.to_excel(writer, sheet_name='Meio_Linha_sem_Comp_TRIF');
    valores_resultados_meio_linha_sem_compensacao_AC.to_excel(writer, sheet_name='Meio_Linha_sem_Comp_AC');
    valores_resultados_meio_linha_sem_compensacao_BC.to_excel(writer, sheet_name='Meio_Linha_sem_Comp_BC');
    valores_resultados_meio_linha_sem_compensacao_AB.to_excel(writer, sheet_name='Meio_Linha_sem_Comp_AB');
    valores_resultados_fim_linha_com_compensacao.to_excel(writer, sheet_name='Fim_Linha_com_Comp');
    valores_resultados_meio_linha_com_compensacao_Trif.to_excel(writer, sheet_name='Meio_Linha_com_Comp_TRIF');
    valores_resultados_meio_linha_com_compensacao_AC.to_excel(writer, sheet_name='Meio_Linha_com_Comp_AC');
    valores_resultados_meio_linha_com_compensacao_BC.to_excel(writer, sheet_name='Meio_Linha_com_Comp_BC');
    valores_resultados_meio_linha_com_compensacao_AB.to_excel(writer, sheet_name='Meio_Linha_com_Comp_AB');

localizacao_Trif = 'DadosTrifasicoT2F.xlsx'
with pd.ExcelWriter(localizacao_Trif) as writer:
    valores_resultados_simulacao_fim_linha.to_excel(writer, sheet_name='Fim_Linha_sem_comp');
    valores_resultados_simulacao_meio_linha_Trif.to_excel(writer, sheet_name='Meio_Linha_TRIF_sem_comp');
    valores_resultados_simulacao_meio_linha_AC.to_excel(writer, sheet_name='Meio_Linha_AC_sem_comp');
    valores_resultados_simulacao_meio_linha_BC.to_excel(writer, sheet_name='Meio_Linha_BC_sem_comp');
    valores_resultados_simulacao_meio_linha_AB.to_excel(writer, sheet_name='Meio_Linha_AB_sem_comp');
    valores_resultados_meio_linha_sem_curto.to_excel(writer, sheet_name='Meio_Linha_Sem_Curto_sem_comp');
    valores_resultados_simulacao_fim_linha_com_comp.to_excel(writer, sheet_name='Fim_Linha_com_comp');
    valores_resultados_simulacao_meio_linha_Trif_com_comp.to_excel(writer, sheet_name='Meio_Linha_TRIF_com_comp');
    valores_resultados_simulacao_meio_linha_AC_com_comp.to_excel(writer, sheet_name='Meio_Linha_AC_com_comp');
    valores_resultados_simulacao_meio_linha_BC_com_comp.to_excel(writer, sheet_name='Meio_Linha_BC_com_comp');
    valores_resultados_simulacao_meio_linha_AB_com_comp.to_excel(writer, sheet_name='Meio_Linha_AB_com_comp');
    valores_resultados_meio_linha_sem_curto_com_comp.to_excel(writer, sheet_name='Meio_Linha_Sem_Curto_com_comp');
    
print('Dados exportados, para '+ localizacao_T2F);
print('Terminado T2F Calculo e Simulação!');
print('Dados exportados, para '+ localizacao_Trif);
print('Terminado T2F Trifásico!');
