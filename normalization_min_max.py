import pandas as pd
import ast

import matplotlib.pyplot as plt
import seaborn as sns

# Min-max normalization
#pLDDT (from AlphaFold) — typically ranges between 0 and 100.
#RMSF (from ATLAS) — varies more freely (often 0–10, but sometimes much higher).
#B-factors (from SoftDis) — can be very large (even 200 in your data).
'''
# Load the CSV file
df = pd.read_csv('data/merged_disorder_scores.csv')

# Columns to normalize
score_columns = ['plddt_alphafold', 'rmsf_atlas', 'bfactors_softdis']

# Safely convert string representations of lists into real Python lists
for col in score_columns:
    df[col] = df[col].apply(ast.literal_eval)

# Expand rows: each residue gets its own row
expanded_rows = []

for _, row in df.iterrows():
    pdb_id = row['pdb_id']
    plddt = row['plddt_alphafold']
    rmsf = row['rmsf_atlas']
    bfac = row['bfactors_softdis']

    length = min(len(plddt), len(rmsf), len(bfac))  # protect against unequal lengths
    for i in range(length):
        expanded_rows.append({
            'pdb_id': pdb_id,
            'residue_index': i + 1,
            'plddt': plddt[i],
            'rmst': rmsf[i],
            'bfactor': bfac[i]
        })

df_long = pd.DataFrame(expanded_rows)
'''
df_long = pd.read_csv("outputs/ATLAS_merged_disorder_scores_long.csv")


# Normalize per pdb_id
def min_max_normalize(group):
    for col in ['plddt', 'rmsf', 'bfactor']:
        min_val = group[col].min()
        max_val = group[col].max()
        if min_val != max_val:
            group[col + '_norm'] = (group[col] - min_val) / (max_val - min_val)
        else:
            group[col + '_norm'] = 0.0
    return group

df_norm = df_long.groupby('pdb_id', group_keys=False).apply(min_max_normalize)

# Save result
df_norm.to_csv('outputs/ATLAS_min_max_normalized_residue_scores.csv', index=False)
print("✅ Normalization complete. Output saved to 'min_max_normalized_residue_scores.csv'.")

df_raw = df_long.copy()

import seaborn as sns
import matplotlib.pyplot as plt


def plot_distributions_all():
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    metrics = ['plddt', 'rmsf', 'bfactor']

    for i, metric in enumerate(metrics):
        #sns.kdeplot(df_raw[metric], ax=axes[i], label='Raw', linewidth=2, linestyle='--')
        sns.kdeplot(df_norm[metric + '_norm'], ax=axes[i], label='Normalized', linewidth=2)
        axes[i].set_title(f'{metric} Distribution (All Proteins)')
        axes[i].set_xlabel('Score')
        axes[i].set_ylabel('Density')
        axes[i].legend()

    plt.tight_layout()
    plt.savefig('plots/ATLAS_distribution_all_proteins_min_max.png')
    plt.show()

plot_distributions_all()

