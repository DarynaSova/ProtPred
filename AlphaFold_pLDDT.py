import pandas as pd

# Load the pLDDT rawdata from ATLAS
df = pd.read_csv("rawdata/plddt_results.csv")

# Group by pdb_id and collect pLDDT values into a list
plddt_by_pdb = df.groupby('pdb_id')['pLDDT'].apply(list).to_dict()

import csv

with open('data/AlphaFold_pLDDT.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['pdb_id', 'pLDDT'])
    for key, value in plddt_by_pdb.items():
        writer.writerow([key, value])
