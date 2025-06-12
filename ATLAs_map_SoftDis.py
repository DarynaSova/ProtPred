import pandas as pd
from datasets import load_from_disk

# load ATLAS files and SoftDis
SoftDis_dataset = load_from_disk("rawdata/SoftDis/clusters_arrow")
SoftDis_data = SoftDis_dataset.to_pandas()
ATLAS_data = pd.read_csv("rawdata/2023_03_09_ATLAS_info.tsv", sep="\t")


# Join on PDB ID
merged_df = pd.merge(SoftDis_data, ATLAS_data, left_on="id", right_on="PDB", how="inner")
columns_to_keep = [
    "id", "UniProt"
]

filtered_df = merged_df[columns_to_keep]
# Save merged results
filtered_df.to_csv("data/ATLAS_SoftDis.csv", index=False)


# Count unique proteins in each table
num_proteins_csv = SoftDis_data["id"].nunique()
num_proteins_tsv = ATLAS_data["PDB"].nunique()
num_proteins_merged = merged_df["id"].nunique()

# Print the results
print(f"Number of unique proteins in SoftDis: {num_proteins_csv}")
print(f"Number of unique proteins in ATLAS: {num_proteins_tsv}")
print(f"Number of unique proteins in merged table: {num_proteins_merged}")

#Number of unique proteins in SoftDis: 26752
#Number of unique proteins in atlas_results: 1938
#Number of unique proteins in merged table: 406


pdb_ids = merged_df["id"]
sequences = merged_df["sequence_x"]

with open("mapping/pdb_codes.txt", "w") as f:
    for pdb_id, sequence in zip(pdb_ids, sequences):
        f.write(pdb_id + "  " + sequence + "\n")

with open("mapping/pdb_codes.fasta", "w") as f:
    for pdb_id, sequence in zip(pdb_ids, sequences):
        f.write(f">{pdb_id}\n{sequence}\n")

with open("rawdata/pdb_codes_for_ATLAS.txt", "w") as f:
    for pdb_id in pdb_ids:
        f.write(pdb_id + "\n")
