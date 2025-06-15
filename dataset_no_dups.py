import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from itertools import product
from itertools import combinations
import os
import re

# This is to remove the instances where there duplicates in the Database 

df = pd.read_csv(r"final_final_filtered_dataset.csv") 
'''
# Count the frequency of each metric in the 'sources' column
metric_counts = df['sources'].value_counts()
print("Category counts:")
print(metric_counts)
'''
duplicate_values = df[df.duplicated(subset=['values'], keep=False)]
print(f"Found {len(duplicate_values)} rows with duplicated 'values'.")
#duplicate_values.to_csv("duplicate_values_only.csv", index=False)

#check for dups 
# Standardize column names (optional if inconsistent)
df.columns = df.columns.str.strip().str.lower()

# Rename for easier handling (you can skip this if already lowercased)
base_cols = ['uniprot_id', 'sources']
# Combine SP_BEG and SP_END into a tuple for pairwise comparison
df['sp_region'] = list(zip(df['sp_beg'], df['sp_end']))
raw_cols = ['pdb_id', 'values', 'sp_region']  # use 'sp_region' directly

# Create output folder
output_dir = "full_duplicates_by_type"
os.makedirs(output_dir, exist_ok=True)

def sanitize_filename(name):
    return re.sub(r'[^\w\d-]', '_', name.lower())

#to avoid overlap
assigned_indices = set()
results = []

# Loop over all non-empty combinations of comparison columns
# Check combinations in descending specificity
for r in range(len(raw_cols), 0, -1):
    for combo in combinations(raw_cols, r):
        label = f"Same {', '.join(base_cols + list(combo))}"
        use_cols = base_cols + list(combo)
        dups = df[df.duplicated(subset=use_cols, keep=False)].copy()
        dups = dups[~dups.index.isin(assigned_indices)]
        if not dups.empty:
            assigned_indices.update(dups.index)
            #add in csv separately
            filename = sanitize_filename(label) + ".csv"
            dups[['pdb_id', 'uniprot_id', 'sources', 'sp_beg', 'sp_end', 'values']].to_csv(os.path.join(output_dir, filename), index=False)
            results.append({
                'Duplicate Type': label,
                'Count': len(dups),
            })

# Display results
summary_df = pd.DataFrame(results)
print(summary_df)



'''
# Filter only 'bfactors' source rows
bfactors_df = df[df['sources'] == 'bfactors'].copy()
print(f"Total rows from 'bfactors': {len(bfactors_df)}")
other_sources_df = df[df['sources'] != 'bfactors'].copy()

# Group by uniprot_id, SP_BEG, SP_END and count how many rows are in each group
duplicate_groups = bfactors_df.groupby(['uniprot_id', 'SP_BEG', 'SP_END']).size().reset_index(name='count')

# Filter only groups with count > 1 (duplicates)
duplicates_only = duplicate_groups[duplicate_groups['count'] > 1]
print(f"Number of duplicate regions \n(same uniprot_id, SP_BEG, SP_END): {len(duplicates_only)}")

# Get full duplicate rows
dup_keys = bfactors_df.duplicated(subset=['uniprot_id', 'SP_BEG', 'SP_END'], keep=False)
duplicate_rows = bfactors_df[dup_keys]
print(f"Number of duplicate rows (all copies): {len(duplicate_rows)}")

# Keep the row with the highest B-factor 'value' for each duplicate group
best_disorder = bfactors_df.sort_values(by='values', ascending=False).drop_duplicates(
    subset=['uniprot_id', 'SP_BEG', 'SP_END']
)

# Logging
print(f"Number of rows after deduplication: {len(best_disorder)}")
print(f"Number of rows removed: {len(bfactors_df) - len(best_disorder)}")

# Combine deduplicated bfactors with all other sources
final_df = pd.concat([best_disorder, other_sources_df], ignore_index=True)

# Save final deduplicated dataset
final_df.to_csv("final_final_nodups.csv", index=False)
print("Saved deduplicated dataset to 'final_final_nodups.csv'")

#final counts all
print(final_df['sources'].value_counts())

'''