import pandas as pd
import ast
import numpy as np

df = pd.read_csv("data/ATLAS_plddt_rmsf_bfactors.csv", quotechar='"', escapechar='\\', encoding='utf-8')
df["pdb_id"] = df["pdb_id"].str.split("_").str[0]

# Clean and convert list-like strings to actual lists for each column
def safe_list_eval(s):
    if isinstance(s, str) and s.strip().startswith("["):
        try:
            # Replace literal nan strings with "nan" so ast can parse it
            s_clean = s.replace('nan', 'np.nan')  # This won't work inside ast.literal_eval directly
            s_clean = s.replace('nan', 'None')    # Better: replace with None
            return ast.literal_eval(s_clean)
        except Exception:
            return np.nan
    return np.nan

# Apply to columns containing list strings
df['bfactors'] = df['bfactors'].apply(safe_list_eval)
df['rmsf_values'] = df['rmsf_values'].apply(safe_list_eval)
df['plddt_values'] = df['plddt_values'].apply(safe_list_eval)
df_filtered = df[df.apply(lambda row: len(set([
    len(row['plddt_values']),
    len(row['rmsf_values']),
    len(row['bfactors'])
])) == 1, axis=1)]

print(len(df), len(df_filtered))
# Recreate the long-format DataFrame with the cleaned dataset
long_df_clean = pd.DataFrame()

for _, row in df_filtered.iterrows():
    length = len(row['plddt_values'])
    temp_df = pd.DataFrame({
        'pdb_id': [row['pdb_id']] * length,
        'residue': list(range(1, length + 1)),
        'plddt': row['plddt_values'],
        'rmsf': row['rmsf_values'],
        'bfactor': row['bfactors']
    })
    long_df_clean = pd.concat([long_df_clean, temp_df], ignore_index=True)

# from 516 proteins only 77 have the same length on Alphafold, fragment needs to be extracted --> sequence needed
long_df_clean.to_csv("outputs/ATLAS_merged_disorder_scores_long.csv")

