import pandas as pd

# Load CSV and TSV
tsv_df1 = pd.read_csv("rawdata/DisProt/DisProt_fixed_new.tsv", sep="\t")
tsv_df2 = pd.read_csv("rawdata/2023_03_09_ATLAS_info.tsv", sep="\t")


# Join on PDB ID
merged_df = pd.merge(tsv_df1, tsv_df2, left_on="acc", right_on="UniProt", how="inner")
columns_to_keep = [
    'acc', "PDB"
]

filtered_df = merged_df[columns_to_keep]
# Save merged results
filtered_df.to_csv("data/DisProt_ATLAS.csv", index=False)


# Count unique proteins in each table
num_proteins_csv = tsv_df1["acc"].nunique()
num_proteins_tsv = tsv_df2["UniProt"].nunique()
num_proteins_merged = merged_df["acc"].nunique()

# Print the results
print(f"Number of unique proteins in CSV: {num_proteins_csv}")
print(f"Number of unique proteins in TSV: {num_proteins_tsv}")
print(f"Number of unique proteins in merged table: {num_proteins_merged}")

#Number of unique proteins in CSV: 3113
#Number of unique proteins in TSV: 1543
#Number of unique proteins in merged table: 118