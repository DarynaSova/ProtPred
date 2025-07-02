import pandas as pd
import ast

# 1. Load the dataset
df = pd.read_csv('dataset/mmseq_final_final_dataset.csv')

# 2. Parse the 'values' column from its string form into a Python list
#    (adjust the column name if yours is different)
df['values_list'] = df['values'].apply(ast.literal_eval)

# 3. Compute residue count per row
df['residue_count'] = df['values_list'].apply(len)

# 4. Sum residue counts by metric (e.g. 'sources')
residue_counts = (
    df
    .groupby('sources')['residue_count']
    .sum()
    .reset_index(name='total_residues')
)

print(residue_counts)

'''
            sources  total_residues
0          bfactors         8657962  no, 767725
1  disprot_disorder         1508454  ok
2            gscore           48147  44066
3             plddt         2257040  ok
4              rmsf          211400  184281
'''