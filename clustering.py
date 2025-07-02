import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.impute import SimpleImputer
import matplotlib.pyplot as plt

# ————— 1) Read in uniprot_id + residue + your three z-score cols —————
df = pd.read_csv(
    'dataset/residue_level_wide_format_norm.csv',
    usecols=['uniprot_id', 'residue_index',
             'bfactors_zscore', 'plddt_zscore', 'rmsf_zscore', 'gscore_zscore']
)

# ————— 2) Make a MultiIndex (uniprot_id, residue) —————
df.set_index(['uniprot_id', 'residue_index'], inplace=True)

# Now df.index is a MultiIndex; df.shape might be something like (N_residues, 3)
# ————— 3) (Optional) Drop any all-NaN columns —————
df = df.dropna(axis=1, how='all')
# 2) Now drop any residue (row) that has a NaN in any of the three metrics:
df_complete = df.dropna(axis=0, how='any')

print(f"Kept {len(df_complete)} residues out of {len(df)} with full data")



# ————— 4) Impute and run PCA as before —————
imputer = SimpleImputer(strategy='mean')
X = imputer.fit_transform(df_complete.values)

pca = PCA(n_components=2)
pcs = pca.fit_transform(X)

# ————— 5) Plot PC1 vs PC2 —————
plt.figure(figsize=(8,8))
plt.scatter(pcs[:, 0], pcs[:, 1], alpha=0.6, s=20)
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.title('PCA Projection by each residue')

plt.tight_layout()
plt.show()
