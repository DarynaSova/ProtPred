import pandas as pd
from datasets import load_from_disk
import numpy as np
import ast
from Bio import SeqIO



def extract_unique_id_positions(file_path):
    """
    Extracts lines with unique IDs from the second column and their positions
    from the 7th and 8th columns of a TSV file.

    Parameters:
        file_path (str): Path to the input TSV file.

    Returns:
        pandas.DataFrame: DataFrame with columns ['Unique_ID', 'Position_Start', 'Position_End']
                          for rows with unique IDs.
    """
    df = pd.read_csv(file_path, sep='\t', header=None)

    # Identify unique values in the second column (index 1)
    unique_ids = df[1].drop_duplicates(keep=False)

    # Filter the DataFrame
    filtered_df = df[df[1].isin(unique_ids)]

    # Extract the relevant columns
    result_df = filtered_df[[1, 6, 7]]
    result_df.columns = ['Unique_ID', 'Position_Start', 'Position_End']

    return result_df

#mmseq_rmsf = extract_unique_id_positions("dataset/uniprot_rmsf_mmseqs.tsv")
#mmseq_bfactors = extract_unique_id_positions("dataset/uniprot_bfactor_mmseqs.tsv")
#mmseq_trizod = extract_unique_id_positions("dataset/uniprot_trizod_mmseqs.tsv")

import pandas as pd
import ast  # safely parses strings of lists into Python list objects


def load_rmsf_to_dict(csv_path):
    """
    Reads a CSV file and returns a dictionary mapping pdb_id to rmsf_values.

    Parameters:
        csv_path (str): Path to the CSV file.

    Returns:
        dict: Dictionary where keys are pdb_ids and values are lists of rmsf_values.
    """
    df = pd.read_csv(csv_path)

    # Convert string representations of lists to actual lists
    rmsf_dict = {
        row['pdb_id'].lower(): ast.literal_eval(row['rmsf_values'])
        for _, row in df.iterrows()
    }
    return rmsf_dict

#rmsf_values = load_rmsf_to_dict("data/ATLAS_plddt_rmsf_bfactors.csv")

def load_bfactors_to_dict(csv_path):
    """
    Reads a CSV file and returns a dictionary mapping pdb_id to rmsf_values.

    Parameters:
        csv_path (str): Path to the CSV file.

    Returns:
        dict: Dictionary where keys are pdb_ids and values are lists of rmsf_values.
    """
    df = pd.read_csv(csv_path)

    # Convert string representations of lists to actual lists
    rmsf_dict = {
        row['pdb_id'].lower(): ast.literal_eval(row['values'])
        for _, row in df.iterrows()
    }
    return rmsf_dict

#bfactors_values = load_bfactors_to_dict("dataset/bfactors_filtered.csv")
#gscores_values = load_bfactors_to_dict("dataset/gscore.csv")



def map_mmseqs_to_values(mmseq_df, value_dict, source_label):
    """
    Maps mmseq alignment segments to structural values from a dictionary.

    Parameters:
        mmseq_df (pd.DataFrame): DataFrame with columns ['Unique_ID', 'Position_Start', 'Position_End']
                                 where 'Unique_ID' is in the form 'uniprot_pdbid'.
        value_dict (dict): Dictionary mapping lowercased pdb_ids to a list of values.
        source_label (str): One of 'rmsf', 'bfactors', 'gscore'.

    Returns:
        list[dict]: Mapped result rows.
    """
    rows = []
    for _, row in mmseq_df.iterrows():
        full_id = row['Unique_ID']
        if '_' not in full_id:
            continue  # skip malformed IDs

        parts = full_id.split('_')

        uniprot_id = parts[0]
        if len(parts) == 2:
            pdb_id = parts[1].lower()
        else:
            pdb_id = parts[1].lower() + "_" + parts[2].lower()
        sp_beg = row['Position_Start']
        sp_end = row['Position_End']
        values = value_dict.get(pdb_id, []) if pdb_id in value_dict else []
        if values:  # Skip empty value ranges
            rows.append({
                'pdb_id': pdb_id,
                'uniprot_id': uniprot_id,
                'sources': source_label,
                'SP_BEG': sp_beg,
                'SP_END': sp_end,
                'values': values
            })
    return rows


# Use previously defined mmseq and value dictionaries (assumed already loaded in context)
#rmsf_data = map_mmseqs_to_values(mmseq_rmsf, rmsf_values, 'rmsf')
#bfactor_data = map_mmseqs_to_values(mmseq_bfactors, bfactors_values, 'bfactors')
#gscore_data = map_mmseqs_to_values(mmseq_trizod, gscores_values, 'gscore')


# Combine and convert to DataFrame
#final_df = pd.DataFrame(rmsf_data + bfactor_data + gscore_data)

# Save to CSV
#output_path = "dataset/mmseq_dataset_rmsf_bfactors_gscores.csv"
#final_df.to_csv(output_path, index=False)


def compile_Disprot(input_combined, output_path):
    # 1) Read in your base dataframe
    mapped_df = pd.read_csv(input_combined)

    # 2) Build a dict of full lengths from your FASTA
    idmapping_path = 'data/idmapping_2025_06_20.fasta'
    full_length = {}
    for record in SeqIO.parse(idmapping_path, "fasta"):
        full_length[record.id.split("|")[1]] = len(record.seq)

    # 3) Load DisProt annotations & filter to 'disorder'
    file_path = 'rawdata/DisProt/DisProt_fixed_new.tsv'
    df = pd.read_csv(file_path, sep='\t')
    disorder_df = df[df['term_name'] == 'disorder'].copy()
    disorder_df['values'] = disorder_df['region_sequence'].apply(
        lambda seq: [True] * len(seq) if pd.notnull(seq) else []
    )

    # 4) Build output rows
    common_columns = ['pdb_id', 'uniprot_id', 'sources', 'SP_BEG', 'SP_END', 'values']
    output_rows = []
    for acc, group in disorder_df.groupby('acc'):
        if acc in full_length:
            length = full_length[acc]
            if length is None:
                continue

            mask = np.zeros(length, dtype=bool)
            for _, row in group.iterrows():
                start = int(row['start']) - 1
                end = int(row['end'])
                mask[start:end] = True


            output_rows.append({
                'pdb_id': None,
                'uniprot_id': acc,
                'sources': 'disprot_disorder',
                'SP_BEG': 1,
                'SP_END': length,
                'values': mask.tolist()
            })

        else:
            output_rows.append({
                'pdb_id': None,
                'uniprot_id': acc,
                'sources': 'disprot_disorder',
                'SP_BEG': row['start'],
                'SP_END': row['end'],
                'values': row['values']
            })



    # 5) Create DataFrame with enforced column order
    disprot_rows = pd.DataFrame(output_rows, columns=common_columns)
    #print(disprot_rows.values)

    # 7) Concatenate, sort, and save
    combined_df = pd.concat([mapped_df, disprot_rows], ignore_index=True)
    combined_df = combined_df.sort_values(by='uniprot_id').reset_index(drop=True)
    combined_df.to_csv(output_path, index=False)
    disprot_rows.to_csv("dataset/disprot.csv", index=False)
    combined_df.to_csv(output_path, index=False)


compile_Disprot("dataset/mmseq_filtered_dataset_rmsf_bfactors_gscores.csv", "dataset/mmseq_final_dataset.csv")

def add_plddt_to_dataset(
    final_csv: str,
    plddt_csv: str,
    output_csv: str,
    id_col: str = 'uniprot_id',
    plddt_col: str = 'plddt_values',
    source_label: str = 'plddt',
    sp_beg: int = 1
) -> pd.DataFrame:
    """
    Reads an existing dataset (final_csv) and a PLDDT file (plddt_csv),
    appends PLDDT rows, saves to output_csv, and returns the combined DataFrame.

    Parameters
    ----------
    final_csv : str
        Path to the original mmseq_final_dataset CSV.
    plddt_csv : str
        Path to the uniprot_plddt_from_fasta CSV.
    output_csv : str
        Path where the combined CSV will be written.
    id_col : str
        Name of the Uniprot ID column in both files.
    plddt_col : str
        Name of the column containing the list of PLDDT values.
    source_label : str
        Value to put in the 'sources' column for new PLDDT rows.
    sp_beg : int
        Starting residue index for PLDDT sequences (usually 1).

    Returns
    -------
    pd.DataFrame
        The concatenated DataFrame with PLDDT rows appended.
    """
    # 1) Load the existing dataset
    df_final = pd.read_csv(final_csv)

    # 2) Load PLDDT data
    df_plddt = pd.read_csv(plddt_csv, usecols=[id_col, plddt_col])

    # 3) Build new rows
    new_rows = []
    for _, row in df_plddt.iterrows():
        uid = row[id_col]
        raw = row[plddt_col]

        # parse string list into Python list
        if isinstance(raw, str):
            try:
                vals = ast.literal_eval(raw)
            except (ValueError, SyntaxError):
                vals = [float(x) for x in raw.strip('[]').split(',') if x]
        else:
            vals = raw if isinstance(raw, list) else list(raw)

        n = len(vals)
        new_rows.append({
            **{c: '' for c in df_final.columns if c not in (id_col, 'sources', 'SP_BEG', 'SP_END', 'values')},
            id_col: uid,
            'sources': source_label,
            'SP_BEG': sp_beg,
            'SP_END': n,
            'values': vals
        })

    # 4) Append and save
    df_new = pd.DataFrame(new_rows, columns=df_final.columns)
    df_out = pd.concat([df_final, df_new], ignore_index=True)
    df_out = df_out.sort_values(by=id_col).reset_index(drop=True)
    df_out.to_csv(output_csv, index=False)
    print(f"Saved combined dataset with {len(df_new)} new PLDDT rows to '{output_csv}'")

    return df_out


add_plddt_to_dataset(final_csv='dataset/mmseq_final_dataset.csv',
                    plddt_csv='dataset/uniprot_plddt_from_fasta.csv',
                    output_csv='dataset/mmseq_final_final_dataset.csv'
)
