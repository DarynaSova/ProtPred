import pandas as pd
import matplotlib.pyplot as plt

# 1) Load everything
df = pd.read_csv('dataset/binary_residue_data_with_disprot.csv')

# 2) Treat every column except 'uniprot_id' and 'residue' as a disorder‐flag metric
metric_cols = [c for c in df.columns if c not in ('uniprot_id', 'residue')]

# 3) Classify
def classify_row(r):
    vals = r[metric_cols].values
    if all(vals == 1):
        return 'always disordered'
    elif all(vals == 0):
        return 'never disordered'
    else:
        return 'partially disordered'

df['classification'] = df.apply(classify_row, axis=1)

# 4) Quick counts & bar‐plot
counts = df['classification'].value_counts().rename_axis('class').reset_index(name='count')
print(counts)

plt.figure(figsize=(6,4))
plt.bar(counts['class'], counts['count'])
plt.xticks(rotation=45)
plt.ylabel('Number of Residues')
plt.title('Residue Disorder Classification (including DisProt flag)')
plt.tight_layout()
plt.show()
