import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# This is to remove the instances where in the bfactor metric there are instances with 
# the same uniprot id and residue start-end but different pdb ids and values. 
# Likely caused by different experimentation on the same region.
# Values kept based on the higher bfactors, as in higher disorder 

# Load the CSV file
df = pd.read_csv(r"final_final_filtered_dataset.csv") 

# Count the frequency of each metric in the 'sources' column
metric_counts = df['sources'].value_counts()
print("Category counts:")
print(metric_counts)

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

