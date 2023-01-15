# Projet_long_PU

This program allows to calculate the number of contact between type of PU on the provided protein dataset

## Packages

Biopython, Numpy, Pandas and Os are required for running this programme

## How to set the inputs

### PDB folder

Put in this folder the pdb files of the proteins whose contacts you want to calculate. These proteins must have 2 chains, "A" and "B".


### OUT folder

Put in this folder the files with the PU calculated for the proteins. For each protein there must be 2 PU files, one for each chain.

The format of these files must be similar to those given as examples.

For the names, the PU files must imperatively be called Out_prot_A or Out_prot_B where prot is the identifier of the protein and is the same as that of the associated pdb file.

## Run the program

You can run this program with batch by typing the commande: python main.py when you are in the PU_contact folder

The program can take a long time to run, the longest being to calculate the inter and intra chain contacts for each protein.
To allow monitoring of the program's status, it will display the protein on which it is calculating contacts and for which chains

## Output files

This program will return a csv file with a table with the contacts for each type of PU. The separators of this file are commas.
