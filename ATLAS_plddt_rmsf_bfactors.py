import os
import pandas as pd


base_dir = "rawdata/atlas_results"
all_data = []

for pdb_id in os.listdir(base_dir):
    pdb_dir = os.path.join(base_dir, pdb_id)
    if not os.path.isdir(pdb_dir):
        continue

    try:
        # Define paths
        plddt_file = os.path.join(pdb_dir, f"{pdb_id}_pLDDT.tsv")
        rmsf_file = os.path.join(pdb_dir, f"{pdb_id}_RMSF.tsv")
        bfactors_file = os.path.join(pdb_dir, f"{pdb_id}_Bfactor.tsv")

        # Read files
        plddt_df = pd.read_csv(plddt_file, sep="\t")
        rmsf_df = pd.read_csv(rmsf_file, sep="\t")
        bfactors_df = pd.read_csv(bfactors_file, sep="\t")

        # Extract numeric values
        plddt_values = plddt_df.iloc[:, 1].astype(float).tolist()
        rmsf_values = rmsf_df.iloc[:, 1:].astype(float).mean(axis=1).tolist()
        bfactors_values = bfactors_df.iloc[:, 0].astype(float).tolist()


        # === Compare lengths
        len_plddt = len(plddt_values)
        len_rmsf = len(rmsf_values)
        len_bfactors = len(bfactors_values)

        consistent = (len_plddt == len_rmsf == len_bfactors)

        if not consistent:
            print(f"Length mismatch for {pdb_id}: pLDDT={len_plddt}, RMSF={len_rmsf}, Bfactors={len_bfactors}")

        # Only include proteins that have bfactors
        if bfactors_df is not None and consistent:
            all_data.append({
                "pdb_id": pdb_id,
                "plddt_values": plddt_values,
                "rmsf_values": rmsf_values,
                "bfactors": bfactors_values,
                "length": len_rmsf
            })
        else:
            print(f"No bfactors found for {pdb_id}")

    except Exception as e:
        print(f"Error processing {pdb_id}: {e}")

# === 3. Save combined rawdata to CSV ===
df = pd.DataFrame(all_data)
df.to_csv("data/ATLAS_plddt_rmsf_bfactors.csv", index=False)