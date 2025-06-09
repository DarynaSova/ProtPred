import pandas as pd
from datasets import load_from_disk
import numpy as np
import ast

def compile_softdis(output_path):
    column_names = ['PDB', 'CHAIN', 'SP_PRIMARY', 'RES_BEG', 'RES_END',
                    'PDB_BEG', 'PDB_END', 'SP_BEG', 'SP_END']

    pdb_uniprot_df = pd.read_csv('mapping/pdb_chain_uniprot.csv', names=column_names, dtype=str)
    pdb_uniprot_df['pdb_id'] = pdb_uniprot_df['PDB'].str.lower() + '_' + pdb_uniprot_df['CHAIN'].str.lower()

    bfactor_dataset = load_from_disk("rawdata/SoftDis/clusters_arrow")
    softdis = bfactor_dataset.to_pandas()
    softdis_clean = softdis[["id", "bfactors"]]
    softdis_clean = softdis_clean.copy()
    softdis_clean['pdb_id'] = softdis_clean['id'].str.lower()
    softdis_clean['bfactors'] = softdis_clean['bfactors'].apply(lambda x: x.tolist() if isinstance(x, np.ndarray) else x)

    print(softdis_clean['bfactors'].head(10))
    print(softdis_clean['bfactors'].apply(type).value_counts())


    softdis_merged = pd.merge(softdis_clean, pdb_uniprot_df, on='pdb_id', how='left')
    softdis_rows = softdis_merged[['pdb_id', 'SP_PRIMARY', 'SP_BEG', 'SP_END']]
    softdis_rows = softdis_rows.copy()
    softdis_rows.rename(columns={'SP_PRIMARY': 'uniprot_id'}, inplace=True)
    softdis_rows['values'] = softdis_merged['bfactors']
    softdis_rows['sources'] = 'bfactors'

    # Reorder columns with 'sources' in 3rd position
    cols = list(softdis_rows.columns)
    cols.remove('sources')
    cols.insert(2, 'sources')
    combined_df = softdis_rows[cols]
    combined_df = combined_df.sort_values(by='pdb_id').reset_index(drop=True)
    combined_df.to_csv(output_path, index=False)


#output_combined_path = "dataset/bfactors.csv"
#compile_softdis(output_combined_path)



def compile_ATLAS_rmsf(output_path):
    column_names = ['PDB', 'CHAIN', 'SP_PRIMARY', 'RES_BEG', 'RES_END',
                    'PDB_BEG', 'PDB_END', 'SP_BEG', 'SP_END']

    # Load the cleaned UniProt-PDB mapping file (skip malformed header row if necessary)
    pdb_uniprot_df = pd.read_csv('mapping/pdb_chain_uniprot.csv', names=column_names, dtype=str)

    # Drop the first row if it's still the header repeated as data
    if pdb_uniprot_df.iloc[0]['PDB'].strip().lower() == 'pdb':
        pdb_uniprot_df = pdb_uniprot_df.drop(index=0).reset_index(drop=True)

    # Construct consistent lowercase pdb_id (e.g., '1ab1_a')
    pdb_uniprot_df['pdb_id'] = pdb_uniprot_df['PDB'].str.lower() + '_' + pdb_uniprot_df['CHAIN'].str.lower()

    # Load the ATLAS data and normalize pdb_id
    atlas_df = pd.read_csv('data/ATLAS_plddt_rmsf_bfactors.csv')
    atlas_df['pdb_id'] = atlas_df['pdb_id'].str.lower()


    # Merge both datasets on pdb_id
    merged_df = pd.merge(atlas_df, pdb_uniprot_df, on='pdb_id', how='left')

    # Select and rename the required columns
    result_df = merged_df[['pdb_id', 'SP_PRIMARY', 'SP_BEG', 'SP_END', 'rmsf_values']]
    result_df = result_df.copy()
    result_df.rename(columns={'SP_PRIMARY': 'uniprot_id'}, inplace=True)
    result_df.rename(columns={'rmsf_values': 'values'}, inplace=True)
    result_df['sources'] = 'RMSF'

    cols = list(result_df.columns)
    cols.remove('sources')
    cols.insert(2, 'sources')
    combined_df = result_df[cols]
    combined_df = combined_df.sort_values(by='pdb_id').reset_index(drop=True)
    combined_df.to_csv(output_path, index=False)

#output_path = "dataset/rmsf.csv"
#compile_ATLAS_rmsf(output_path)

def compile_ATLAS_plddt(output_path):
    column_names = ['PDB', 'CHAIN', 'SP_PRIMARY', 'RES_BEG', 'RES_END',
                    'PDB_BEG', 'PDB_END', 'SP_BEG', 'SP_END']

    # Load the cleaned UniProt-PDB mapping file (skip malformed header row if necessary)
    pdb_uniprot_df = pd.read_csv('mapping/pdb_chain_uniprot.csv', names=column_names, dtype=str)

    # Drop the first row if it's still the header repeated as data
    if pdb_uniprot_df.iloc[0]['PDB'].strip().lower() == 'pdb':
        pdb_uniprot_df = pdb_uniprot_df.drop(index=0).reset_index(drop=True)

    # Construct consistent lowercase pdb_id (e.g., '1ab1_a')
    pdb_uniprot_df['pdb_id'] = pdb_uniprot_df['PDB'].str.lower() + '_' + pdb_uniprot_df['CHAIN'].str.lower()

    # Load the ATLAS data and normalize pdb_id
    atlas_df = pd.read_csv('data/ATLAS_plddt_rmsf_bfactors.csv')
    atlas_df['pdb_id'] = atlas_df['pdb_id'].str.lower()


    # Merge both datasets on pdb_id
    merged_df = pd.merge(atlas_df, pdb_uniprot_df, on='pdb_id', how='left')

    # Select and rename the required columns
    result_df = merged_df[['pdb_id', 'SP_PRIMARY', 'SP_BEG', 'SP_END', 'plddt_values']]
    result_df = result_df.copy()
    result_df.rename(columns={'SP_PRIMARY': 'uniprot_id'}, inplace=True)
    result_df.rename(columns={'plddt_values': 'values'}, inplace=True)
    result_df['sources'] = 'plddt'

    cols = list(result_df.columns)
    cols.remove('sources')
    cols.insert(2, 'sources')
    combined_df = result_df[cols]
    combined_df = combined_df.sort_values(by='pdb_id').reset_index(drop=True)
    combined_df.to_csv(output_path, index=False)

output_path = "dataset/plddt.csv"
compile_ATLAS_plddt(output_path)

import pandas as pd

def combine_and_sort_csv(rmsf_path, plddt_path, bfactors_path, output_path):
    """
    Combine RMSF and PLDDT CSVs and sort by uniprot_id.

    Parameters:
        rmsf_path (str): Path to the RMSF CSV file
        plddt_path (str): Path to the PLDDT CSV file
        bfactors_path (str): Path to the bfactors CSV file
        output_path (str): Path to save the combined and sorted CSV
    """
    # Load input files
    rmsf_df = pd.read_csv(rmsf_path)
    plddt_df = pd.read_csv(plddt_path)
    bfactors_df = pd.read_csv(bfactors_path)

    # Combine the two datasets
    combined_df = pd.concat([rmsf_df, plddt_df, bfactors_df], ignore_index=True)
    cols = list(combined_df.columns)
    cols.remove('sources')
    cols.insert(2, 'sources')
    combined_df = combined_df[cols]
    # Sort by uniprot_id
    combined_df_sorted = combined_df.sort_values(by='uniprot_id').reset_index(drop=True)

    combined_df_sorted['values'] = combined_df_sorted['values'].astype(str).str.replace('\n', ' ', regex=True)
    combined_df_sorted.to_csv(output_path, index=False)


#if __name__ == "__main__":
combine_and_sort_csv(
        rmsf_path="dataset/rmsf.csv",
        plddt_path="dataset/plddt.csv",
        bfactors_path="dataset/bfactors.csv",
       output_path="dataset/combined_rmsf_plddt_bfactors.csv"
    )


def compile_Disprot_gscore(input_combined, output_path):
    mapped_df = pd.read_csv(input_combined)

    # Load the TSV file
    file_path = 'rawdata/DisProt/DisProt_fixed_new.tsv'
    df = pd.read_csv(file_path, sep='\t')

    # Filter the dataframe to keep only the specified columns
    filtered_df = df[['acc', 'start', 'end', 'term_name', 'region_sequence']]

    # Filter rows where term_name is "disorder"
    disorder_df = df[df['term_name'] == 'disorder'].copy()

    # Create the new column 'values' with a list of True values based on the length of 'region_sequence'
    disorder_df['values'] = disorder_df['region_sequence'].apply(lambda seq: [True] * len(seq) if pd.notnull(seq) else [])

    common_columns = ['pdb_id', 'uniprot_id', 'sources', 'SP_BEG', 'SP_END', 'values']


    disprot_rows = pd.DataFrame({
        'pdb_id': None,
        'uniprot_id': disorder_df['acc'],
        'sources': 'disprot_disorder',
        'SP_BEG': disorder_df['start'],
        'SP_END': disorder_df['end'],
        'values': disorder_df['values']
    })[common_columns]

    mapped_df = mapped_df[common_columns]

    # Ensure column order matches
    disprot_rows = disprot_rows[mapped_df.columns]

    # Append and save
    combined_df = pd.concat([mapped_df, disprot_rows], ignore_index=True)
    combined_df = combined_df.sort_values(by='uniprot_id').reset_index(drop=True)
    mapped_df = combined_df
    gscore_df = pd.read_csv("data/merged_trizod_bmrb_output.csv")
    # Load the PDB-UniProt mapping file
    column_names = ['PDB', 'CHAIN', 'SP_PRIMARY', 'RES_BEG', 'RES_END', 'PDB_BEG', 'PDB_END', 'SP_BEG', 'SP_END']
    mapping_df = pd.read_csv('mapping/pdb_chain_uniprot.csv', names=column_names, skiprows=1, dtype=str)

    # Drop any rows that are malformed
    mapping_df = mapping_df[mapping_df['PDB'].str.len() == 4]

    # Construct lowercase pdb_id to match the gscore file
    mapping_df['pdb_id'] = mapping_df['PDB'].str.lower()
    gscore_df['pdb_id'] = gscore_df['pdb_id'].str.lower()

    merged_gscore = pd.merge(gscore_df, mapping_df, on='pdb_id', how='left')

    # Prepare the new rows
    gscore_rows = pd.DataFrame({
        'pdb_id': merged_gscore['pdb_id'],
        'uniprot_id': merged_gscore['SP_PRIMARY'],
        'sources': 'gscore',
        'SP_BEG': merged_gscore['SP_BEG'],
        'SP_END': merged_gscore['SP_END'],
        'values': merged_gscore['gscore']
    })

    gscore_rows = gscore_rows[mapped_df.columns]

    # Append and save
    combined_df = pd.concat([mapped_df, gscore_rows], ignore_index=True)
    combined_df = combined_df.sort_values(by='uniprot_id').reset_index(drop=True)
    combined_df = combined_df.drop_duplicates(subset=['pdb_id', 'sources', 'SP_BEG'])

    combined_df.to_csv(output_path, index=False)

compile_Disprot_gscore("dataset/combined_rmsf_plddt_bfactors.csv", "dataset/final_final_dataset.csv")

