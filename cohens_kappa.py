
import pandas as pd
from sklearn.metrics import cohen_kappa_score
import ast
import numpy as np

# Path to your input CSV
INPUT_CSV = "dataset/residue_level_wide_format_zscore_norm.csv"
df = pd.read_csv(INPUT_CSV)

# Threshold definitions


thresholds = {
    'bfactors_zscore': (2, 'greater'),
    'plddt_zscore': (0.125, 'less'),       # 👈 reversed logic
    'rmsf_zscore': (11.477, 'greater'),
    'gscore_zscore': (2.640, 'greater')
}
'''

# best threshold
thresholds = {
    'bfactors_zscore': (1.9889, 'greater'),
    'plddt_zscore': (-0.4137, 'less'),       # 👈 reversed logic
    'rmsf_zscore': (1.5164, 'greater'),
    'gscore_zscore': (0.3974    , 'greater')
}
'''

zscore_cols = [col for col in df.columns if col.endswith('_zscore') and col in thresholds]


def threshold_to_binary(series, threshold, direction='greater'):
    return series.apply(lambda x: np.nan if pd.isna(x) else int(x > threshold) if direction == 'greater' else int(x < threshold))

binary_columns = {}

for col in thresholds:
    threshold, direction = thresholds[col]  # ✅ unpack tuple
    binary_columns[col + '_bin'] = threshold_to_binary(df[col], threshold, direction)

# Include disprot_disorder if present
if 'disprot_disorder' in df.columns:
    binary_columns['disprot_disorder'] = df['disprot_disorder'].astype('float')

# --- Manual Cohen's kappa function ---
def compute_cohen_kappa_manual(x, y):
    tp = np.sum((x == 1) & (y == 1))
    tn = np.sum((x == 0) & (y == 0))
    fp = np.sum((x == 0) & (y == 1))
    fn = np.sum((x == 1) & (y == 0))
    total = tp + tn + fp + fn
    po = (tp + tn) / total if total != 0 else 0
    px_yes = (tp + fn) / total if total != 0 else 0
    py_yes = (tp + fp) / total if total != 0 else 0
    px_no = (tn + fp) / total if total != 0 else 0
    py_no = (tn + fn) / total if total != 0 else 0
    pe = px_yes * py_yes + px_no * py_no
    return (po - pe) / (1 - pe) if (1 - pe) != 0 else 0


# --- Pairwise Cohen's kappa computation ---
results = {}
binary_keys = list(binary_columns.keys())

for i, col1 in enumerate(binary_keys):
    for col2 in binary_keys[i + 1:]:
        bin1 = binary_columns[col1]
        bin2 = binary_columns[col2]
        pair_df = pd.DataFrame({col1: bin1, col2: bin2}).dropna()
        if not pair_df.empty:
            kappa = compute_cohen_kappa_manual(pair_df[col1].values, pair_df[col2].values)
            results[(col1, col2)] = kappa
        else:
            results[(col1, col2)] = None

print("check value counts")
for i, col in enumerate(binary_keys):
    print(f"{col}:{binary_columns[col].value_counts().to_dict()}")

# --- Save results ---
kappa_df = pd.DataFrame.from_dict(results, orient='index', columns=['Cohen Kappa (NaN-aware)'])
kappa_df.to_csv("data/pairwise_cohen_kappa_with_disprot.csv")

binary_df = pd.DataFrame(binary_columns)

# Optionally add identifying columns like uniprot_id and residue_index if needed
binary_df['uniprot_id'] = df['uniprot_id']
binary_df['residue_index'] = df['residue_index']

# --- After building binary_columns dict ---

# Count zeros and ones for each binary metric
counts = {}
for col, series in binary_columns.items():
    # drop NaN so they don't get counted
    non_na = series.dropna().astype(int)
    counts[col] = {
        'count_0': int((non_na == 0).sum()),
        'count_1': int((non_na == 1).sum())
    }

# Turn into a DataFrame
counts_df = pd.DataFrame.from_dict(counts, orient='index')

# Print to console
print("\nZero/One counts per metric:")
print(counts_df)

# Optionally save to CSV
counts_df.to_csv("data/binary_counts_per_metric.csv")

'''
Zero/One counts per metric original threshold:
                     count_0  count_1
bfactors_zscore_bin   667508   100217
plddt_zscore_bin     1443760   813280
rmsf_zscore_bin       184253       28
gscore_zscore_bin      41634     2432
disprot_disorder     1261095   247359


Zero/One counts per metric best threshold:
                     count_0  count_1
bfactors_zscore_bin   667342   100383
plddt_zscore_bin     1610697   646343
rmsf_zscore_bin       173625    10656
gscore_zscore_bin      37933     6133
disprot_disorder     1261095   247359
'''




# --- Save to CSV ---
binary_df.to_csv("dataset/binary_residue_data_with_disprot_first_threshold.csv", index=False)

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Raw Cohen's kappa data
data = {
    ('bfactors', 'plddt'): 0.5730931359066979,
    ('bfactors', 'rmsf'): 0.007396354757963077,
    ('bfactors', 'gscore'): 0.29032944627408025,
    ('bfactors', 'disprot'): 0.5924837340066369,
    ('plddt', 'rmsf'): 0.0014089122266112096,
    ('plddt', 'gscore'): 0.24268012903170696,
    ('plddt', 'disprot'): 0.2786257077885735,
    ('rmsf', 'gscore'): 0.0,
    ('rmsf', 'disprot'): 0.0018377338331817489,
    ('gscore', 'disprot'): 0.3290903108262953,
}

# Extract unique labels
labels = sorted({lab for pair in data.keys() for lab in pair})

# Create an empty matrix and fill it
matrix = np.eye(len(labels))
label_index = {label: idx for idx, label in enumerate(labels)}

for (a, b), value in data.items():
    i, j = label_index[a], label_index[b]
    matrix[i, j] = value
    matrix[j, i] = value

# Create a heatmap
fig, ax = plt.subplots()
cax = ax.imshow(matrix)

# Configure ticks and labels
ax.set_xticks(np.arange(len(labels)))
ax.set_yticks(np.arange(len(labels)))
ax.set_xticklabels(labels, rotation=45, ha='right')
ax.set_yticklabels(labels)

# Annotate with values
for i in range(len(labels)):
    for j in range(len(labels)):
        ax.text(j, i, f"{matrix[i, j]:.2f}", ha="center", va="center")

plt.title("Pairwise Cohen's Kappa Heatmap")
plt.tight_layout()
plt.show()



'''
G scores threshold: 0.65

'''
'''
Bfacotrs: 83
'''

'''
PLDDT: 80
Regions with pLDDT > 90 are expected to be modelled to high accuracy. These should be suitable for any application that benefits from high accuracy (e.g. characterising binding sites).
Regions with pLDDT between 70 and 90 are expected to be modelled well (a generally good backbone prediction).
Regions with pLDDT between 50 and 70 are low confidence and should be treated with caution.
'''

'''
RMSF: 2
We defined residues with mean square deviations in the NMR ensemble of greater than 2Å as disordered; this definition was supported by visual inspection of the ensembles.
Wang RY, Han Y, Krassovsky K, Sheffler W, Tyka M, Baker D. Modeling disordered regions in proteins using Rosetta. PLoS One. 2011;6(7):e22060. doi: 10.1371/journal.pone.0022060. Epub 2011 Jul 29. PMID: 21829444; PMCID: PMC3146542.
'''