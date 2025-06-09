import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load the CSV file
df = pd.read_csv(r"C:\Users\Dasha\ProtPred\mapped_all_sources_with_disprot.csv") 

# Parse the 'values' column (assuming it's a string of comma-separated numbers)
def parse_values(val):
    if isinstance(val, str):
        val = val.replace('[', '').replace(']', '').strip()
        floats = []
        for x in val.split(','):
            x = x.strip()
            try:
                floats.append(float(x))
            except ValueError:
                continue  # Skip non-numeric values
        return np.array(floats)
    elif isinstance(val, (list, np.ndarray)):
        return np.array(val)
    else:
        return np.array([])

df['parsed_values'] = df['values'].apply(parse_values)

# Normalize values using Min-Max normalization
def min_max_normalize(arr):
    if len(arr) == 0:
        return arr
    min_val = np.min(arr)
    max_val = np.max(arr)
    if max_val == min_val:
        return np.zeros_like(arr)
    return (arr - min_val) / (max_val - min_val)

df['normalized_values'] = df['parsed_values'].apply(min_max_normalize)

# Example: Visualize normalized values for the first 3 entries
plt.figure(figsize=(10, 6))
for idx, row in df.head(3).iterrows():
    plt.plot(row['normalized_values'], label=f"Entry {idx} ({row.get('type', 'unknown')})")
plt.xlabel('Index')
plt.ylabel('Normalized Value')
plt.title('Normalized Values for First 3 Entries')
plt.legend()
plt.tight_layout()
plt.show()

# Save normalized data to a new CSV (flatten normalized values as strings)
df['normalized_values_str'] = df['normalized_values'].apply(lambda arr: ','.join(map(str, arr)))
df.to_csv('normalized_output.csv', index=False)

print("Normalization complete. Output saved to 'normalized_output.csv'.")