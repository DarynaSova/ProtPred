import pandas as pd
import ast
from datasets import load_from_disk
import csv

'''
# === 1. Load B-factors from Arrow Dataset ===
bfactor_dataset = load_from_disk("rawdata/SoftDis/clusters_arrow")
df = bfactor_dataset.to_pandas()

df["bfactors"] = df["bfactors"].apply(lambda x: x.tolist() if hasattr(x, "tolist") else x)
df.to_csv("data/SoftDis_bfactors_clean.csv", index=False)


'''

from datasets import load_from_disk
import pandas as pd
import json

# Load your Arrow dataset
dataset = load_from_disk("rawdata/SoftDis/clusters_arrow")

# Convert to pandas DataFrame
df = dataset.to_pandas()

# Ensure list-like fields (e.g., bfactors) are serialized properly
# json.dumps will turn [89.49, 99.93, 46.59] → "[89.49, 99.93, 46.59]"
df["bfactors"] = df["bfactors"].apply(json.dumps)

# Save to a clean JSON file (one record per line, or all as a list)
df.to_json("data/bfactors_clean.json", orient="records", indent=2)
