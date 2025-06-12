import pandas as pd
import ast

df = pd.read_csv("dataset/final_final_dataset.csv")

# get values from different sources in different dataframe and transfer them to lists
df_pldddt = df[df['sources'] == 'PLDDT']
df_pldddt['values'] = df_pldddt['values'].apply(ast.literal_eval)

df_bfactors = df[df['sources'] == 'bfactors']
df_bfactors['values'] = df_bfactors['values'].apply(ast.literal_eval)

df_RMSF = df[df['sources'] == 'RMSF']
df_RMSF['values'] = df_RMSF['values'].apply(ast.literal_eval)

df_gscore = df[df['sources'] == 'gscore']
df_gscore['values'] = df_gscore['values'].apply(ast.literal_eval)

# input the normalized threshold for each metrices
th_plddt =
th_bfactors =
th_RMSF =
th_gscore =


