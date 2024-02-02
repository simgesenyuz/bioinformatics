# Bioinformatics

This repo contains notebooks related to bioinformatics/computational biology that I have written for various tasks.

## <ins>Mapping Uniprot Residues to PDB/mmCIF Residues:</ins> <br />
&nbsp; In this task, we have the following dataset: [Docking Results](Residue_Mapping/Docking_Results.txt). This dataset contains PDB chains and the task is to get all residues of the PDB and map UniProt residues to PDB residues. <br /> 
&nbsp; Related codes: [Residue Mapping](Residue_Mapping/map_uniprot_residues_to_pdb_residues.ipynb)<br />  <br /> 
&nbsp; <ins> Skills Used:</ins> <br /> 
&nbsp; &nbsp; data cleaning with pandas DataFrames <br /> 
&nbsp; &nbsp; using APIs to collect data in JSON format without saving any file anywhere. Specific APIs used: [PDBe API](https://www.ebi.ac.uk/pdbe/api/) and [PDBe Graph API](https://www.ebi.ac.uk/pdbe/graph-api/pdbe_doc/)

## <ins>PDB Extract Ligand/Drug Binding Pockets v1:</ins> <br />
&nbsp; In this task, we identify binding pockets within protein structures using PDB codes.

&nbsp; Features:  <br /> 
&nbsp; PDB Code Input: Users can input PDB codes directly to analyze specific protein structures (given PDB id and chain id). <br /> 
&nbsp; Automated Pocket Extraction: Identifies potential ligand/drug binding sites automatically. <br /> 
&nbsp; Related codes: [Pocket Extraction](Scripts_Drug_Repurposing/pdb_extract_pockets_v1.ipynb)<br />  <br /> 
&nbsp; <ins> Skills Used:</ins> <br /> 
&nbsp; &nbsp; using REST APIs to post and collect data in JSON format Specific API used: [ProteinPlus API](https://proteins.plus/help/dogsite_rest)
