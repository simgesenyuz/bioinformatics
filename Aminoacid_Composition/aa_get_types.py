import sys, os
import pandas as pd
from ast import literal_eval
from collections import Counter


def get_res_type(res):
    """
    This function groups residues based on hydrophobicity etc.
    :param res: list of residues
    :return: resdiue counts based on properties, dict type
    """
    residue_categories = {
        'ALA': 'Hydrophobic', 'VAL': 'Hydrophobic', 'LEU': 'Hydrophobic', 'ILE': 'Hydrophobic',
        'MET': 'Hydrophobic', 'PHE': 'Hydrophobic', 'TRP': 'Hydrophobic', 'PRO': 'Hydrophobic',
        'GLN': 'Polar', 'ASN': 'Polar', 'SER': 'Polar',
        'THR': 'Polar', 'TYR': 'Hydrophobic', 'CYS': 'Polar', 'GLY': 'Nonpolar',
        'LYS': 'Positively Charged', 'ARG': 'Positively Charged', 'HIS': 'Positively Charged',
        'ASP': 'Negatively Charged', 'GLU': 'Negatively Charged'
    }

    res_categories = [residue_categories.get(residue, 'Unknown') for residue in res]
    res_counts = dict(Counter(res_categories))
    return res_counts


os.chdir('../../aa_prop')
df = pd.DataFrame()

# Save all related files from a directory
for filename in os.listdir():
    if filename.startswith('aa_comp'):
        temp_df = pd.read_csv(filename,  sep='\t', index_col=None)
        df = df._append(temp_df, ignore_index=True)

# If the column has a value of a str, represent the value as a list, else save an empty list
df.loc[:, 'res_chain_1_contacting'] = df['res_chain_1_contacting'].apply(lambda x: literal_eval(x) if type(x) == str else list())
df.loc[:, 'res_chain_1_nearby'] = df['res_chain_1_nearby'].apply(lambda x: literal_eval(x) if type(x) == str else list())
df.loc[:, 'res_chain_2_contacting'] = df['res_chain_2_contacting'].apply(lambda x: literal_eval(x) if type(x) == str else list())
df.loc[:, 'res_chain_2_nearby'] = df['res_chain_2_nearby'].apply(lambda x: literal_eval(x)if type(x) == str else list())

# Run get_res_type() function on data to get residue types
df.loc[:, 'type_res_chain_1_contacting'] = df['res_chain_1_contacting'].apply(get_res_type)
df.loc[:, 'type_res_chain_1_nearby'] = df['res_chain_1_nearby'].apply(get_res_type)
df.loc[:, 'type_res_chain_2_contacting'] = df['res_chain_2_contacting'].apply(get_res_type)
df.loc[:, 'type_res_chain_2_nearby'] = df['res_chain_2_nearby'].apply(get_res_type)

# Save to file
df.to_csv('aa_with_properties.txt', sep='\t')
