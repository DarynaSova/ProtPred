import pandas as pd
import matplotlib.pyplot as plt
from upsetplot import UpSet, from_memberships

# Load the data
df = pd.read_csv("c:/Users/Dasha/ProtPred/residue_level_wide_format.csv")

# Define disorder using thresholds
df['disordered_plddt'] = df['PLDDT'] < 0.2
df['disordered_rmsf'] = df['RMSF'] > 0.8
df['disordered_bfactor'] = df['bfactors'] > 0.8

# Create a list of memberships: which methods call each residue disordered
memberships = []
for _, row in df.iterrows():
    methods = []
    if bool(row['disordered_plddt']):
        methods.append("PLDDT")
    if bool(row['disordered_rmsf']):
        methods.append("RMSF")
    if bool(row['disordered_bfactor']):
        methods.append("BFACTOR")
    if methods:
        memberships.append(methods)

print(memberships[:10])  # Debug: check structure

# Check for empty memberships
if not memberships:
    raise ValueError("No disordered residues found. Check your thresholds or input data.")

# Build the UpSet data
upset_data = from_memberships(memberships)

plt.figure(figsize=(10, 6))
upset = UpSet(upset_data, subset_size='count', show_counts=True)
upset.plot()
plt.suptitle("Overlap of Disorder Definitions\n", fontsize=14)
plt.tight_layout()
plt.show()