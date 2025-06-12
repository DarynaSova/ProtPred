
import pandas as pd
import ast



# Path to your input CSV
INPUT_CSV = "dataset/final_final_nodups.csv"
# Path for output
OUTPUT_CSV = "dataset/final_binary_values.csv"

# Threshold definitions
THRESHOLDS = {
    'bfactors': 83,
    'plddt': 80,
    'RMSF': 2,
    'G scores': 0.65
}


def to_binary(source, parsed_values):
    """
    Convert a list of floats into a list of booleans according to:
      - If source in THRESHOLDS:  v > threshold
      - If source == 'disprot_disorder': bool(v)  (pass-through)
      - Otherwise:               None (skip)
    """
    if source in THRESHOLDS:
        thresh = THRESHOLDS[source]
        return [v > thresh for v in parsed_values]
    elif source == 'disprot_disorder':
        # assume values are 0/1 already
        return [bool(v) for v in parsed_values]
    else:
        return None


def main():
    # 1) Load
    df = pd.read_csv(INPUT_CSV)

    # 2) Parse the 'values' column (string) into Python lists
    df['parsed_values'] = df['values'].apply(ast.literal_eval)

    # 3) Compute binary_values
    df['binary_values'] = df.apply(
        lambda row: to_binary(row['sources'], row['parsed_values']),
        axis=1
    )

    # 4) Drop rows we couldn’t handle (if any)
    df_clean = df.dropna(subset=['SP_BEG', 'binary_values'])
    df_clean = df_clean.drop(columns=['parsed_values'])
    # 5) Save result
    df_clean.to_csv(OUTPUT_CSV, index=False)
    print(f"Wrote {len(df_clean)} rows with binary_values to {OUTPUT_CSV}")


if __name__ == "__main__":
    main()




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