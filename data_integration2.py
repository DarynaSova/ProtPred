# Re-import necessary libraries after code execution environment reset
import pandas as pd
import ast

# Reload the file
file_path = 'outputs/merged_disorder_scores.csv'
df = pd.read_csv(file_path)

# Clean and convert list-like strings to actual lists for each column
df['plddt_alphafold'] = df['plddt_alphafold'].apply(lambda x: ast.literal_eval(x))
df['rmsf_atlas'] = df['rmsf_atlas'].apply(lambda x: ast.literal_eval(x))
df['bfactors_softdis'] = df['bfactors_softdis'].apply(lambda x: list(map(float, x.strip('[]').split())))

# Filter out rows with mismatched lengths
df_filtered = df[df.apply(lambda row: len(set([
    len(row['plddt_alphafold']),
    len(row['rmsf_atlas']),
    len(row['bfactors_softdis'])
])) == 1, axis=1)]

print(print(len(df)), len(df_filtered))
# Recreate the long-format DataFrame with the cleaned dataset
long_df_clean = pd.DataFrame()

for _, row in df_filtered.iterrows():
    length = len(row['plddt_alphafold'])
    temp_df = pd.DataFrame({
        'pdb_id': [row['pdb_id']] * length,
        'residue': list(range(1, length + 1)),
        'plddt': row['plddt_alphafold'],
        'rmsf': row['rmsf_atlas'],
        'bfactor': row['bfactors_softdis']
    })
    long_df_clean = pd.concat([long_df_clean, temp_df], ignore_index=True)

# from 516 proteins only 77 have the same length on Alphafold, fragment needs to be extracted --> sequence needed
long_df_clean.to_csv("outputs/merged_disorder_scores_long.csv")