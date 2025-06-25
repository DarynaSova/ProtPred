import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from sklearn.metrics import precision_score, recall_score, f1_score
from sklearn.metrics import cohen_kappa_score
import seaborn as sns


# --- Load dataset ---
df = pd.read_csv(r"residue_level_wide_format_zscore_norm.csv")

# --- Manual Cohen's Kappa Function ---
def compute_cohen_kappa_manual(x, y):
    tp = np.sum((x == 1) & (y == 1))
    tn = np.sum((x == 0) & (y == 0))
    fp = np.sum((x == 1) & (y == 0))  
    fn = np.sum((x == 0) & (y == 1))  
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

# --- Updated function to optimize threshold and track metrics ---
def optimize_threshold_against_disprot(series, disprot, direction='greater', steps=100):
    valid = ~series.isna() & ~disprot.isna()
    series = series[valid]
    disprot = disprot[valid]

    # Clip threshold scan to the z-score range of approximately [-3, 3]
    clipped_min = max(series.min(), -3)
    clipped_max = min(series.max(), 3)
    thresholds = np.linspace(clipped_min, clipped_max, steps)
    
    best_kappa = -1
    best_threshold = None
    best_metrics = {}
    metrics_per_threshold = []

    for t in thresholds:
        if direction == 'greater':
            binary = series > t
        else:
            binary = series < t

        y_pred = binary.astype(int).values
        y_true = disprot.astype(int).values

        kappa = compute_cohen_kappa_manual(y_pred, y_true)
        precision = precision_score(y_true, y_pred, zero_division=0)
        recall = recall_score(y_true, y_pred, zero_division=0)
        f1 = f1_score(y_true, y_pred, zero_division=0)

        metrics_per_threshold.append({
            'threshold': t,
            'kappa': kappa,
            'precision': precision,
            'recall': recall,
            'f1_score': f1
        })

        if kappa > best_kappa:
            best_kappa = kappa
            best_threshold = t
            best_metrics = {
                'precision': precision,
                'recall': recall,
                'f1_score': f1
            }

    return best_threshold, best_kappa, best_metrics, metrics_per_threshold

# --- Predictors and threshold directions ---
predictors = {
    'bfactors_zscore': 'greater',
    'plddt_zscore': 'less',
    'rmsf_zscore': 'greater',
    'gscore_zscore': 'greater'
}

# --- Prepare output directory for plots ---
os.makedirs("plots", exist_ok=True)

# --- Run optimization and plot results ---
results = []

for col, direction in predictors.items():
    if col in df.columns and 'disprot_disorder' in df.columns:
        threshold, kappa, metrics, metric_trace = optimize_threshold_against_disprot(
            df[col], df['disprot_disorder'], direction=direction, steps=100
        )

        # Save main summary result
        results.append({
            'metric': col,
            'direction': direction,
            'best_threshold': threshold,
            'kappa': kappa,
            'precision': metrics['precision'],
            'recall': metrics['recall'],
            'f1_score': metrics['f1_score']
        })

        # Plot metrics vs threshold
        trace_df = pd.DataFrame(metric_trace)
        plt.figure(figsize=(10, 6))
        plt.plot(trace_df['threshold'], trace_df['kappa'], label='Kappa')
        plt.plot(trace_df['threshold'], trace_df['f1_score'], label='F1 Score')
        plt.plot(trace_df['threshold'], trace_df['precision'], label='Precision', linestyle='--')
        plt.plot(trace_df['threshold'], trace_df['recall'], label='Recall', linestyle='--')
        plt.axvline(x=threshold, color='gray', linestyle=':', label='Best Threshold')
        plt.xlabel('Threshold')
        plt.ylabel('Metric Value')
        plt.title(f'Threshold Sensitivity for {col}')
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(f'plots/{col}_threshold_sensitivity.png')
        plt.close()

# --- Save results ---
results_df = pd.DataFrame(results)
results_df.to_csv("optimized_thresholds_vs_disprot.csv", index=False)
#print(results_df[['metric', 'best_threshold']])

# --- Output summary ---
print(results_df)
print("Saved results to 'optimized_thresholds_vs_disprot.csv'")
print("Plots saved in 'plots/' directory.")


'''
#Best thresholds Cohen's Kappa corellation
# Binarize predictors with best thresholds
binary_calls = pd.DataFrame(index=df.index)

for _, row in results_df.iterrows():
    metric = row['metric']
    threshold = row['best_threshold']
    direction = row['direction']

    print(metric, threshold)

    if metric in df.columns:
        if direction == 'greater':
            binary_calls[metric] = (df[metric] > threshold).astype(int)
        else:
            binary_calls[metric] = (df[metric] < threshold).astype(int)

# Calculate pairwise Cohen's kappa
metrics = list(binary_calls.columns)
n = len(metrics)
kappa_matrix = np.zeros((n, n))

for i in range(n):
    for j in range(n):
        valid = binary_calls[[metrics[i], metrics[j]]].notna().all(axis=1)
        y1 = binary_calls.loc[valid, metrics[i]]
        y2 = binary_calls.loc[valid, metrics[j]]
        if len(y1) > 0:
            kappa_matrix[i, j] = cohen_kappa_score(y1, y2)
        else:
            kappa_matrix[i, j] = np.nan

# Plot heatmap
plt.figure(figsize=(8, 6))
sns.heatmap(kappa_matrix, xticklabels=metrics, yticklabels=metrics,
            annot=True, fmt=".2f", cmap="coolwarm", center=0)
plt.title("Pairwise Cohen's Kappa Between Binary Disorder Predictions")
plt.tight_layout()
plt.show()
'''
