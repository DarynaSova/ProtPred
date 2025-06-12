import pandas as pd
import ast
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import zscore

# Load the CSV file
df = pd.read_csv('dataset/final_final_dataset.csv')

df_pldddt = df[df['sources'] == 'PLDDT']
df_pldddt['values'] = df_pldddt['values'].apply(ast.literal_eval)

df_bfactors = df[df['sources'] == 'bfactors']
df_bfactors['values'] = df_bfactors['values'].apply(ast.literal_eval)

df_RMSF = df[df['sources'] == 'RMSF']
df_RMSF['values'] = df_RMSF['values'].apply(ast.literal_eval)

df_gscore = df[df['sources'] == 'gscore']
df_gscore['values'] = df_gscore['values'].apply(ast.literal_eval)

# Expand rows: one row per residue
expanded_rows = []

for single_df in [df_gscore, df_RMSF, df_bfactors, df_pldddt]
    for _, row in single_df.iterrows():
        pdb_id = row['pdb_id']
        values = row['values']

        for i in range(len(values)):
            expanded_rows.append({
                'pdb_id': pdb_id,
                'residue_index': i + 1,
                'values': values[i]
            })
    df_long = pd.DataFrame(expanded_rows)



# Z-score normalization across all proteins (global)
df_long['pLDDT_z'] = zscore(df_long['pLDDT'])
df_long['RMSF_z'] = zscore(df_long['RMSF'])
df_long['Bfactor_z'] = zscore(df_long['Bfactor'])

# Save z-normalized data
df_long.to_csv('outputs/z_normalized_residue_scores.csv', index=False)
print("✅ Z-score normalization complete. Output saved to 'z_normalized_residue_scores.csv'.")

# Compute correlation matrix (global, across all residues and proteins)
correlation_matrix = df_long[['pLDDT_z', 'RMSF_z', 'Bfactor_z']].corr()

# Plot correlation heatmap
plt.figure(figsize=(6, 5))
sns.heatmap(correlation_matrix, annot=True, cmap='coolwarm', fmt=".2f", square=True)
plt.title("Pearson Correlation Heatmap (Z-Normalized Features)")
plt.tight_layout()
plt.savefig('plots/zscore_correlation_heatmap.png')
plt.show()


def plot_distributions_raw_vs_zscore():
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    metric_pairs = [
        ('pLDDT', 'pLDDT_z'),
        ('RMSF', 'RMSF_z'),
        ('Bfactor', 'Bfactor_z')
    ]

    for i, (raw_col, z_col) in enumerate(metric_pairs):
        sns.kdeplot(df_long[raw_col], ax=axes[i], label='Raw', linewidth=2, linestyle='--')
        sns.kdeplot(df_long[z_col], ax=axes[i], label='Z-score Normalized', linewidth=2)

        axes[i].set_title(f'{raw_col} Distribution (Raw vs. Z-score)')
        axes[i].set_xlabel('Value')
        axes[i].set_ylabel('Density')
        axes[i].legend()
        axes[i].grid(True)

        if 'z' in z_col:
            axes[i].axvline(0, color='gray', linestyle=':', linewidth=1)

    plt.tight_layout()
    plt.savefig('plots/raw_vs_zscore_distributions.png')
    plt.show()


# Call the function
plot_distributions_raw_vs_zscore()



