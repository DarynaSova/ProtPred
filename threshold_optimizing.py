import pandas as pd
import numpy as np

# --- Load dataset ---
df = pd.read_csv("dataset/residue_level_wide_format_new_nodups_zscore.csv")

# --- Manual Cohen's Kappa Function ---
def compute_cohen_kappa_manual(x, y):
    tp = np.sum((x == 1) & (y == 1))
    tn = np.sum((x == 0) & (y == 0))
    fp = np.sum((x == 0) & (y == 1))
    fn = np.sum((x == 1) & (y == 0))
    total = tp + tn + fp + fn
    if total == 0:
        return 0
    po = (tp + tn) / total
    px_yes = (tp + fn) / total
    py_yes = (tp + fp) / total
    px_no = (tn + fp) / total
    py_no = (tn + fn) / total
    pe = px_yes * py_yes + px_no * py_no
    return (po - pe) / (1 - pe) if (1 - pe) != 0 else 0

# --- Function to optimize threshold for a column ---
def optimize_threshold_against_disprot(series, disprot, direction='greater', steps=100):
    # Remove rows with missing values
    valid = ~series.isna() & ~disprot.isna()
    series = series[valid]
    disprot = disprot[valid]

    best_kappa = -1
    best_threshold = None

    min_val, max_val = series.min(), series.max()
    thresholds = np.linspace(min_val, max_val, steps)

    for t in thresholds:
        if direction == 'greater':
            binary = series > t
        else:
            binary = series < t

        kappa = compute_cohen_kappa_manual(binary.astype(int).values, disprot.astype(int).values)
        if kappa > best_kappa:
            best_kappa = kappa
            best_threshold = t

    return best_threshold, best_kappa

# --- Predictors and threshold directions ---
predictors = {
    'bfactors_zscore': 'greater',
    'plddt_zscore': 'less',
    'RMSF_zscore': 'greater',
    'gscore_zscore': 'greater'
}

# --- Run optimization ---
results = []

for col, direction in predictors.items():
    if col in df.columns and 'disprot_disorder' in df.columns:
        threshold, kappa = optimize_threshold_against_disprot(
            df[col], df['disprot_disorder'], direction=direction, steps=100
        )
        results.append({'metric': col, 'direction': direction, 'best_threshold': threshold, 'kappa': kappa})

# --- Save or print results ---
results_df = pd.DataFrame(results)
results_df.to_csv("optimized_thresholds_vs_disprot.csv", index=False)

print(results_df)
print("✅ Saved results to 'optimized_thresholds_vs_disprot.csv'")
