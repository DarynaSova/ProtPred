import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


# Load the uploaded CSV file
file_path = "/Users/gaoyiting/Bioinformatik/pp1/create_dataset/result/Softdis_atlas_alphafold_trizod/final_final_nodups.csv"
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
            "RMSF": np.nan,
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


# Z-score normalization
for col in ["bfactors", "RMSF", "plddt", "gscore"]:
    if col in final_df.columns:
        mean = final_df[col].mean(skipna=True)
        std = final_df[col].std(skipna=True)
        if std and not np.isnan(std):
            final_df[col + "_zscore"] = (final_df[col] - mean) / std
        else:
            final_df[col + "_zscore"] = np.nan


final_df.to_csv("/Users/gaoyiting/Bioinformatik/pp1/create_dataset/result/Softdis_atlas_alphafold_trizod/residue_level_wide_format_new_nodups.csv")


corr_matrix = final_df[["plddt_zscore", "bfactors_zscore", "RMSF_zscore", "gscore_zscore"]].corr()


plt.figure(figsize=(8, 6))
sns.heatmap(corr_matrix, annot=True, cmap="coolwarm", fmt=".2f", vmin=-1, vmax=1)
plt.title("Correlation Matrix Between Different Metrics")
plt.tight_layout()
plt.savefig("/Users/gaoyiting/Bioinformatik/pp1/create_dataset/result/correlation_matrix_new.png", dpi=300, bbox_inches='tight')


from itertools import combinations

sources = ["plddt_zscore", "bfactors_zscore", "RMSF_zscore", "gscore_zscore"]
for s1, s2 in combinations(sources, 2):
    n = final_df.dropna(subset=[s1, s2]).shape[0]
    print(f"{s1} ∩ {s2}: {n} residues with both values")