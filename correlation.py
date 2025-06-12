import pandas as pd


df_norm = pd.read_csv('outputs/ATLAS_min_max_normalized_residue_scores.csv')

correlation_matrix = df_norm[['plddt', 'rmsf', 'bfactor',
                              'plddt_norm', 'rmsf_norm', 'bfactor_norm']].corr()

print(correlation_matrix)

import seaborn as sns
import matplotlib.pyplot as plt

plt.figure(figsize=(8, 6))
sns.heatmap(correlation_matrix, annot=True, cmap='coolwarm', fmt=".2f")
plt.title('Pairwise Correlation Between Disorder Scores')
plt.tight_layout()
plt.savefig('plots/ATLAS_correlation_matrix_heatmap.png')
plt.show()
