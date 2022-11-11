import sys
import argparse
import warnings
import Bio.PDB
from Bio import pairwise2
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBIO import PDBIO
from Bio import SeqIO
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.Structure import Structure
from Bio.PDB.Selection import unfold_entities
import pandas as pd

import os

# SOURCE:   https://proteopedia.org/wiki/index.php/Amino_Acids
d3to1 = {'ALA': 'A', 'ASX': 'B', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', \
        'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', \
        'LEU': 'L', 'MET': 'M', 'ASN': 'N', 'PYL': 'O', 'PRO': 'P', \
        'GLN': 'Q', 'ARG': 'R', 'SER': 'S', 'THR': 'T', 'SEC': 'U', \
        'VAL': 'V', 'TRP': 'W', 'UNK': 'X', 'TYR': 'Y', 'GLX': 'Z'}


def contains_sequence(pdb_structure, sequence_1letter):
    """asd
        Args:

        Returns:
    asd
    """
    model_residues = unfold_entities(pdb_structure, 'R')

    model_sequence_1letter = ''.join([d3to1[residue.resname] for residue in model_residues])

    if sequence_1letter in model_sequence_1letter:
        return True
    else:
        return False



def main():

    pdbp = PDBParser(QUIET=True)

    abspath = os.getcwd()
    abspath = abspath.split("GPCR_project")[0]+"GPCR_project"

    pdb_models = abspath + "/data/model_pdb"
    csv_file = abspath + "/data/data.csv"
    data = pd.read_csv(csv_file, dtype="str", sep=',', header=0)

    n = 0
    p = 0
    for row in data.itertuples():
        antigen_sequence = getattr(row, "Antigen_Sequence")
        uniprot_id = getattr(row, "Uniprot_ID")
        pdb_file = pdb_models + "/AF-{}-F1-model_v4.pdb".format(uniprot_id)
        contains_agseq = contains_sequence(pdbp.get_structure("", pdb_file), antigen_sequence)
        if contains_agseq:
            p += 1
        else: 
            n +=1
        print("Uniprot: {}, Contains: {}".format(uniprot_id,contains_agseq))
    print("Positive: {}/{}".format(p,n+p))

if __name__ == '__main__':
    main()