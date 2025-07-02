import matplotlib.pyplot as plt
import numpy as np

# Data
categories = ['B factors', 'RMSF', 'PLDDT', 'G Score', 'Disprot']

# Original thresholds
orig_ord =    [667508, 184253, 1443760,  41634, 1261095]
orig_dis =    [100217,     28,   813280,   2432,   247359]

# Best thresholds
best_ord =    [667342, 173625, 1610697,  37933, 1261095]
best_dis =    [100383,  10656,   646343,   6133,   247359]

x = np.arange(len(categories))
total_width = 0.8
bar_width = total_width / 4

# New color scheme: dark/light blue, dark/light orange
colors = ['#1f77b4',  # Orig Ordered (blue)
          '#AEC7E8',  # Best Ordered (light blue)
          '#ff7f0e',  # Orig Disordered (orange)
          '#ffbb78']  # Best Disordered (light orange)

labels = ['Orig Ordered', 'Best Ordered', 'Orig Disordered', 'Best Disordered']
datasets = [orig_ord, best_ord, orig_dis, best_dis]
offsets = [-1.5, -0.5, 0.5, 1.5]

fig, ax = plt.subplots(figsize=(10, 6))
bars = []
for off, data, label, color in zip(offsets, datasets, labels, colors):
    bars.append(ax.bar(x + off*bar_width, data, bar_width, label=label, color=color))

def autolabel(bars, dx):
    """
    bars: BarContainer from ax.bar()
    dx:   horizontal shift in data coords (positive→right, negative→left)
    """
    ymax = max(orig_ord + best_ord + orig_dis + best_dis)
    for bar in bars:
        h = bar.get_height()
        y = h + 0.01 * ymax + dx
        x = bar.get_x() + bar.get_width() / 2
        ax.text(
            x,
            y,
            f'{h:,}',
            ha='center', va='bottom',
            fontsize=6
        )


turn = False
for bg in bars:
    if turn == False:
        autolabel(bg, dx = 1)
        turn = True
    else:
        autolabel(bg, dx = 2)
        turn = False


ax.set_xticks(x)
ax.set_xticklabels(categories, rotation=45, ha='right')
ax.set_xlabel('Source Metric')
ax.set_ylabel('Residue Count')
ax.set_title('Original vs Best Thresholds (Ordered & Disordered)')
ax.legend(ncol=2, fontsize=9)

plt.tight_layout()
plt.show()



# Data
'''
categories = ['B factors', 'RMSF', 'PLDDT', 'G Score', 'Disprot']
ordered_counts = [667508, 184253, 1443760, 41634, 1261095]
disordered_counts = [100217, 28, 813280, 2432, 247359]



categories = ['B factors', 'RMSF', 'PLDDT', 'G Score', 'Disprot']
ordered_counts = [667342, 173625, 1610697, 37933, 1261095]
disordered_counts = [100383, 10656, 646343, 6133, 247359]



group_spacing = 0.8
x = np.arange(len(categories)) * group_spacing
width = 0.35

fig, ax = plt.subplots(figsize=(10, 6))

# Plot the bars
bars_ord = ax.bar(x - width/2, ordered_counts, width, label='Ordered', color='#1f77b4')
bars_dis = ax.bar(x + width/2, disordered_counts, width, label='Disordered', color='#ff7f0e')

# Add labels on top of each bar
def autolabel(bars):
    for bar in bars:
        height = bar.get_height()
        ax.text(
            bar.get_x() + bar.get_width()/2,  # x-coordinate: center of the bar
            height + 0.01*max(ordered_counts+disordered_counts),  # a little above the top
            f'{height:,}',  # format with commas
            ha='center',
            va='bottom',
            fontsize= 9
        )

autolabel(bars_ord)
autolabel(bars_dis)

# Labels, ticks, legend
ax.set_xlabel('Metric')
ax.set_ylabel('Residue Count')
ax.set_title('Ordered vs Disordered Residue Counts')
ax.set_xticks(x)
ax.set_xticklabels(categories, rotation=45, ha='right')
ax.legend()

plt.tight_layout()
plt.show()

'''

'''

# Data provided by the user
labels = ['bfactors', 'plddt', 'disprot', 'rmsf', 'gscore']
values = [8657962, 2257040, 1508454, 211400, 48147]

# Create bar plot
plt.figure(figsize=(8, 5))
bars = plt.bar(labels, values, color = "blue")
plt.bar(labels, values)
plt.xlabel('Metric')
plt.ylabel('Residue Count')
#plt.title('Residue Counts by Source')
plt.xticks(rotation=45, ha='right')
# Label each bar with its height
for bar in bars:
    height = bar.get_height()
    plt.text(
        bar.get_x() + bar.get_width() / 2,  # x-coordinate
        height,                             # y-coordinate
        f'{height:,}',                      # formatted height
        ha='center',
        va='bottom'
    )

plt.tight_layout()
plt.show()


# Load the dataset
df = pd.read_csv('dataset/mmseq_final_final_dataset.csv')

# Compute total residues per source as sum of list lengths in 'values'
df['residue_count'] = df['values'].apply(lambda x: len(eval(x)) if isinstance(x, str) else len(x))

# Aggregate by source
residue_counts = df.groupby('sources')['residue_count'].sum().reset_index()

print(residue_counts)





def check_value_lengths(csv_path):
    """
    Checks if the length of the 'values' list equals SP_END - SP_BEG + 1.
    Prints rows with mismatched lengths.

    Parameters:
        csv_path (str): Path to the mapped output CSV file.
    """
    df = pd.read_csv(csv_path)

    # Safely convert the 'values' column from string to list
    df['values'] = df['values'].apply(ast.literal_eval)

    for _, row in df.iterrows():
        expected_len = row['SP_END'] - row['SP_BEG'] + 1
        actual_len = len(row['values'])
        if expected_len != actual_len:
            print(f"Length mismatch for pdb_id: {row['pdb_id']}, "
                  f"uniprot_id: {row['uniprot_id']}, source: {row['sources']}, "
                  f"Expected: {expected_len}, Found: {actual_len}")

check_value_lengths("dataset/mmseq_dataset_rmsf_bfactors_gscores.csv")

import pandas as pd
import ast

def filter_valid_rows(csv_path, output_path):
    """
    Filters out rows where the length of 'values' does not match SP_END - SP_BEG + 1,
    and saves the valid rows to a new CSV.

    Parameters:
        csv_path (str): Path to the input CSV file.
        output_path (str): Path to save the filtered CSV file.
    """
    df = pd.read_csv(csv_path)
    df['values'] = df['values'].apply(ast.literal_eval)

    def is_valid(row):
        return len(row['values']) == (row['SP_END'] - row['SP_BEG'] + 1)

    valid_df = df[df.apply(is_valid, axis=1)]
    valid_df.to_csv(output_path, index=False)


#filter_valid_rows("dataset/mmseq_dataset_rmsf_bfactors_gscores.csv", "dataset/mmseq_filtered_dataset_rmsf_bfactors_gscores.csv")

'''

