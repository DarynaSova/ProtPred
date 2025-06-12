import pandas as pd

# Load the uploaded CSV file
file_path = "dataset/final_final_nodups.csv"
df = pd.read_csv(file_path)

import ast
import numpy as np

# Filter out rows without proper SP_BEG and SP_END
df_filtered = df.dropna(subset=["SP_BEG", "SP_END"]).copy()

# Ensure SP_BEG and SP_END are integers
df_filtered["SP_BEG"] = df_filtered["SP_BEG"].astype(int)
df_filtered["SP_END"] = df_filtered["SP_END"].astype(int)

# Convert string representation of lists to actual lists
df_filtered["values"] = df_filtered["values"].apply(ast.literal_eval)

# Initialize an empty list to collect transformed data
records = []

# Iterate through each row to expand it into residue-level data
for _, row in df_filtered.iterrows():
    uniprot_id = row["uniprot_id"]
    source = row["sources"]
    start = row["SP_BEG"]
    values = row["values"]

    for i, val in enumerate(values):
        residue_index = start + i
        record = {
            "uniprot_id": uniprot_id,
            "residue_index": residue_index,
            "bfactors": np.nan,
            "rmsf": np.nan,
            "plddt": np.nan,
            "gscore": np.nan
        }
        record[source] = val  # assign the value to the correct column based on source
        records.append(record)

# Convert the list of dictionaries into a DataFrame
expanded_df = pd.DataFrame(records)

# Group by unique identifiers and aggregate to combine bfactors, rmsf, and plddt for each residue
final_df = expanded_df.groupby(
    ["uniprot_id", "residue_index"], as_index=False
).first()

final_df.to_csv("dataset/residue_level_wide_format_new.csv")