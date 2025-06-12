import pandas as pd
from collections import Counter

# 1) List of file‐paths
file_paths = ['dataset/bfactors.csv', 'dataset/disprot.csv', 'dataset/gscore.csv', 'dataset/plddt.csv', 'dataset/rmsf.csv']

# 2) Count in how many distinct files each ID appears
id_counter = Counter()
for path in file_paths:
    df = pd.read_csv(path, usecols=['uniprot_id'])
    # use set(df['id']) so duplicates *within* a file only count once
    id_counter.update(set(df['uniprot_id']))

# 3) IDs appearing in >=2 files
shared_ids = {i for i, cnt in id_counter.items() if cnt >= 2}

filtered_dfs = []
for path in file_paths:
    df = pd.read_csv(path)
    # keep only rows whose id is in at least one other file
    df = df[df['uniprot_id'].isin(shared_ids)]
    filtered_dfs.append(df)
    # (optional) save out the filtered CSV
    df.to_csv(path.replace('.csv', '_filtered.csv'), index=False)

combined = pd.concat(filtered_dfs, ignore_index=True)
combined.to_csv('dataset/final_final_filtered_dataset.csv', index=False)
