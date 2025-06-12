import pandas as pd
import ast

import matplotlib.pyplot as plt
import seaborn as sns

# Min-max normalization
#pLDDT (from AlphaFold) — typically ranges between 0 and 100.
#RMSF (from ATLAS) — varies more freely (often 0–10, but sometimes much higher).
#B-factors (from SoftDis) — can be very large (even 200 in your data).

# Load the CSV file
df = pd.read_csv('data/ATLAS_plddt_rmsf_bfactors.csv')

# Columns to normalize
score_columns = ['plddt_values', 'rmsf_values', 'bfactors']


# Expand rows: each residue gets its own row
expanded_rows = []

for _, row in df.iterrows():
    pdb_id = row['pdb_id']
    plddt = row['plddt_values'].replace("[", "").replace("]", "").split(", ")
    rmsf = row['rmsf_values'].replace("[", "").replace("]", "").split(", ")
    bfac = row['bfactors'].replace("[", "").replace("]", "").split(", ")

    length = min(len(plddt), len(rmsf), len(bfac))  # protect against unequal lengths

    for i in range(length):
        expanded_rows.append({
            'pdb_id': pdb_id,
            'residue_index': i + 1,
            'pLDDT': float(plddt[i]),
            'RMSF': float(rmsf[i]),
            'Bfactor': float(bfac[i])
        })

df_long = pd.DataFrame(expanded_rows)

# Normalize per pdb_id
def min_max_normalize(group):
    for col in ['pLDDT', 'RMSF', 'Bfactor']:
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
    metrics = ['pLDDT', 'RMSF', 'Bfactor']

    for i, metric in enumerate(metrics):
        sns.kdeplot(df_raw[metric], ax=axes[i], label='Raw', linewidth=2, linestyle='--')
        #sns.kdeplot(df_norm[metric + '_norm'], ax=axes[i], label='Normalized', linewidth=2)
        axes[i].set_title(f'{metric} Distribution (All Proteins)')
        axes[i].set_xlabel('Score')
        axes[i].set_ylabel('Density')
        axes[i].legend()

    plt.tight_layout()
    plt.savefig('plots/ATLAS_distribution_all_proteins_raw.png')
    plt.show()

plot_distributions_all()

