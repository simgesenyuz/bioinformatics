#!/usr/bin/env python
# coding: utf-8

import os, sys
import pandas as pd
import requests
from Bio import SeqIO
import numpy as np
import json


def get_pdb_residues(pdb_id, chain_id):
    """
    :param pdb_id: str
    :param chain_id: str
    :return: pdb_residues, residues of a given pdb_id, chain_id pair
    """
    # get residues of a pdb
    try:
        pdb_id = pdb_id.lower()
        response = requests.get(f'https://www.ebi.ac.uk/pdbe/api/pdb/entry/residue_listing/{pdb_id}/chain/{chain_id}')
        pdb_residues = json.loads(response.content)[pdb_id]['molecules']

        response.close()
        return pdb_residues

    except:
        return None


def get_residue_type(pdb_residues):
    """
    :param pdb_residues: pdb
    :return: residues_dict: dict of residue_name: author_residue_number pairs
    """

    res_list = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 
                'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 
                'THR', 'TRP', 'TYR', 'VAL']

    # Get residue number: residue name pairs as a dictionary
    for entry in pdb_residues:
        for chain in entry['chains']:
            if any(res in [d['residue_name'] for d in chain['residues']] for res in res_list):
                res_numbers = [d['author_residue_number'] for d in chain['residues']]
                res_names = [d['residue_name'] for d in chain['residues']]
                residues_dict = {k: v for k, v in zip(res_numbers, res_names)}
    return residues_dict


def match_res_type(interface_res_list, residues_dict):
    """
    :return: :param interface_res_list: list of residue_numbers (int)
    :param residues_dict: dict of residue_name: author_residue_number pairs
    :return: filtered residue_names
    """
    # Match residue names with a set of residue numbers
    return [residues_dict[res] for res in interface_res_list]


def save_res_type(row, type_interface_res):
    """
    :param row: row of a df, this df contains pdbID (str), chain (str), chain_1_contacting (str)...
    :param type_interface_res: type of interface residues ex. 'chain_1_contacting' (str)
    :return:
    """
    try:
        if len(row[type_interface_res])>0:
            # Extract pdbID and chain_id from the row
            pdbID, chain_id = row['pdbID'], row['chain']

            pdb_residues = get_pdb_residues(pdbID, chain_id) 
            residues_dict = get_residue_type(pdb_residues, chain_id)

            # List of interface residues
            interface_res_list = row[type_interface_res]
    
            # Match the residue types
            res_type_list = match_res_type(interface_res_list, residues_dict)
            
        else:
            res_type_list = []
        
        return res_type_list
    
    except:
        return ''


# Main
# Get filename as sys argument
filename = sys.argv[1]
df = pd.read_csv(filename, sep='\t')

# Use this line if you want to randomly sample from df
df = df.sample(n=5)

# Convert string of contacting/nearby residues to lists
df = df.fillna('')
df.loc[:, 'chain_1_contacting'] = df['chain_1_contacting'].apply(lambda x: [int(res) for res in x.split(', ')] if len(x)>0 else [])
df.loc[:, 'chain_1_nearby'] = df['chain_1_nearby'].apply(lambda x: [int(res) for res in x.split(', ')] if len(x)>0 else [])
df.loc[:, 'chain_2_contacting'] = df['chain_2_contacting'].apply(lambda x: [int(res) for res in x.split(', ')] if len(x)>0 else [])
df.loc[:, 'chain_2_nearby'] = df['chain_2_nearby'].apply(lambda x: [int(res) for res in x.split(', ')] if len(x)>0 else [])

# Save residues
df.loc[:, 'res_chain_1_contacting'] = df.loc[:].apply(lambda row: save_res_type(row, 'chain_1_contacting'), axis=1)
df.loc[:, 'res_chain_1_nearby'] = df.loc[:].apply(lambda row: save_res_type(row, 'chain_1_nearby'), axis=1)
df.loc[:, 'res_chain_2_contacting'] = df.loc[:].apply(lambda row: save_res_type(row, 'chain_2_contacting'), axis=1)
df.loc[:, 'res_chain_2_nearby'] = df.loc[:].apply(lambda row: save_res_type(row, 'chain_2_nearby'), axis=1)

# Save df result
output_filename = filename.split('_')[1].split('.')[0]
df.to_csv(f'aa_comp_{output_filename}.txt', sep='\t', index=None)
