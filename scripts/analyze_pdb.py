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

import numpy as np
import statistics
from Bio.PDB.DSSP import DSSP


import os

# SOURCE:   https://proteopedia.org/wiki/index.php/Amino_Acids
d3to1 = {'ALA': 'A', 'ASX': 'B', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', \
        'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', \
        'LEU': 'L', 'MET': 'M', 'ASN': 'N', 'PYL': 'O', 'PRO': 'P', \
        'GLN': 'Q', 'ARG': 'R', 'SER': 'S', 'THR': 'T', 'SEC': 'U', \
        'VAL': 'V', 'TRP': 'W', 'UNK': 'X', 'TYR': 'Y', 'GLX': 'Z'}

# SS-Scheme 1: H,G,I->H ; E,B->E ; -,T,S->C
dssp8to3 = {'H': 'H', 'G': 'H', 'I': 'H', 'E': 'E', 'B': 'E', \
        'T': 'C', 'S': 'C', '-': 'C'}


def perfect_match(pdb_chain, sequence_1letter: str):
    """ 
        Args:
            pdb_chain: biopython chain object
            sequence_1letter: 1 letter amino acid

        Returns:
            True if input sequence is part of input chain, Flase otherwise.
    
    """
    model_residues = unfold_entities(pdb_chain, 'R')
    model_sequence_1letter = ''.join([d3to1[residue.resname] for residue in model_residues])

    if sequence_1letter in model_sequence_1letter:
        return True
    else:
        return False


def compute_aligned_sequence(pdb_chain, subsequence: str):
    """Align input subsequence to chain sequence.

        Args:
        pdb_chain: biopython chain structure
        subsequence: a 1letter protein sequence

        Returns:
            Returns a subsequence of the pdb_chain, aligned to the input subsequence as well \
                 as a "normalized alignment score" (alignment score divided by length of subsequence) \
                    measuring quality of the alignment (score 1 means subsequence matches part of the input chain)
    """
    model_residues = unfold_entities(pdb_chain, 'R')

    model_sequence = ''.join([d3to1[residue.resname] for residue in model_residues])

    alignments = pairwise2.align.localms(model_sequence, subsequence, 1, -0.5, -1, -0.5)
    # Identical characters are given 1 points,0.5 point is deducted for each non-identical character.
    # 1 points are deducted when opening a gap, and 0.5 points are deducted when extending it
    
    alignment = alignments[0]


    normalized_alignment_score = round(alignment[2]/len(subsequence),2)

    start_pos = alignment[3]
    end_pos = alignment[4]
    aligned_model_sequence = alignment[0][start_pos:end_pos]
    #aligned_subsequence = alignment[1][start_pos:end_pos]

    return aligned_model_sequence, normalized_alignment_score


def get_matched_residue_ids(pdb_chain, subsequence):
    """
        Args:

        Returns:
    
    """

    model_residues = unfold_entities(pdb_chain, 'R')
    model_sequence = ''.join([d3to1[residue.resname] for residue in model_residues])

    if model_sequence.count(subsequence) > 1:
        sys.exit("WARNING: The subsequence appears more than once in the model sequence.")
    if model_sequence.count(subsequence) == 0:
        sys.exit("WARNING: The subsequence does not exist in the model sequence.")

    start = model_sequence.index(subsequence)
    end = start + len(subsequence)
    matched_residue_ids = [residue.get_id() for residue in model_residues[start:end]]
    
    return matched_residue_ids



def compute_relative_plDDT(pdb_chain, subsequence="", threshold=0.25):
    """
        Args:

        Returns:
    
    """

    if subsequence == "":
        model_atoms = unfold_entities(pdb_chain, 'A')
        b_factor_scores = [atom.get_bfactor() for atom in model_atoms]
        relative_plDDT_score = sum([1 for score in b_factor_scores if score <= threshold])/len(b_factor_scores)
    else:
        matched_ids = get_matched_residue_ids(pdb_chain, subsequence)
        subsequence_match = [pdb_chain[residue_id] for residue_id in matched_ids]

        b_factor_scores = [atom.get_bfactor() for residue in subsequence_match for atom in residue.get_atoms()]
        relative_plDDT_score = sum([1 for score in b_factor_scores if score <= threshold])/len(b_factor_scores)
        
    relative_plDDT_score = round(relative_plDDT_score,2)
    return relative_plDDT_score

def compute_average_plDDT(pdb_chain, subsequence=""):
    """
        Args:

        Returns:
    
    """

    if subsequence == "":
        model_atoms = unfold_entities(pdb_chain, 'A')
        b_factor_scores = [atom.get_bfactor() for atom in model_atoms]
        average_plDDT_score = statistics.mean(b_factor_scores)
    else:
        matched_ids = get_matched_residue_ids(pdb_chain, subsequence)
        subsequence_match = [pdb_chain[residue_id] for residue_id in matched_ids]
        b_factor_scores = [atom.get_bfactor() for residue in subsequence_match for atom in residue.get_atoms()]
        average_plDDT_score = statistics.mean(b_factor_scores)
        
    average_plDDT_score = round(average_plDDT_score,2)
    return average_plDDT_score

    """
    else:
        model_residues = unfold_entities(pdb_chain, 'R')
        model_sequence = ''.join([d3to1[residue.resname] for residue in model_residues])

        if model_sequence.count(subsequence) > 1:
            sys.exit("WARNING: The subsequence appears more than once in the model sequence.")
        if model_sequence.count(subsequence) == 0:
            sys.exit("WARNING: The subsequence does not exist in the model sequence.")

        start = model_sequence.index(subsequence)
        end = start + len(subsequence)
        subsequence_match = model_residues[start:end]


        average_plDDT_score = statistics.mean([atom.get_bfactor() for residue in subsequence_match for atom in residue.get_atoms()])
        return average_plDDT_score
    """


def main():


    pdbp = PDBParser(QUIET=True)

    abspath = os.getcwd()
    abspath = abspath.split("GPCR_project")[0]+"GPCR_project"

    pdb_models = abspath + "/data/model_pdb"
    csv_file = abspath + "/data/data.csv"
    csv_out = abspath + "/analysis/data_results.csv"
    data = pd.read_csv(csv_file, dtype="str", sep=',', header=0)

    # initialize new columns
    data["matched_subsequence"] = np.nan
    data["normalized_alignment_score"] = np.nan
    data["average_plDDT_score"] = np.nan
    data["plDDT_below_025"] = np.nan
    data["secondary_structure"] = np.nan
    data["average_relative_ASA"] = np.nan
    data["relative_ASA_below_025"] = np.nan
    data["relative_H-dssp"] = np.nan
    data["relative_E-dssp"] = np.nan
    data["relative_C-dssp"] = np.nan


    for row in data.itertuples():
        index = getattr(row, "Index")
        antigen_sequence = getattr(row, "Antigen_Sequence")
        uniprot_id = getattr(row, "Uniprot_ID")
        pdb_file = pdb_models + "/AF-{}-F1-model_v4.pdb".format(uniprot_id)
        pdb_structure = pdbp.get_structure("", pdb_file)
        model_chain = pdb_structure[0]["A"]

        # Give update
        print("Processing {}/{}".format(index+1, data.shape[0]))

        matched_subsequence, normalized_alignment_score = compute_aligned_sequence(model_chain, antigen_sequence)
        data.at[index,"matched_subsequence"] = matched_subsequence
        data.at[index,"normalized_alignment_score"] = normalized_alignment_score
        data.at[index,"average_plDDT_score"] = compute_average_plDDT(model_chain,matched_subsequence)
        data.at[index,"plDDT_below_025"] = compute_relative_plDDT(model_chain,matched_subsequence,0.25)


        matched_residue_ids = get_matched_residue_ids(model_chain, compute_aligned_sequence(model_chain, antigen_sequence)[0])
        
        # see: https://biopython.org/docs/1.75/api/Bio.PDB.DSSP.html
        dssp = DSSP(pdb_structure[0], pdb_file, dssp='mkdssp', acc_array='Wilke')
        #dssp = DSSP(pdb_structure[0], pdb_file, dssp='mkdssp', acc_array='Sander')

        # analyse secondary structure
        secondary_structure_3letter = ''.join([dssp8to3[dssp[(model_chain.get_id(), residue_id)][2]] for residue_id in matched_residue_ids])
        data.at[index, "secondary_structure"] = secondary_structure_3letter

        data.at[index, "relative_H-dssp"] = round(secondary_structure_3letter.count("H")/len(secondary_structure_3letter),2)
        data.at[index, "relative_E-dssp"] = round(secondary_structure_3letter.count("E")/len(secondary_structure_3letter),2)
        data.at[index, "relative_C-dssp"] = round(secondary_structure_3letter.count("C")/len(secondary_structure_3letter),2)

        # analyse relative ASA
        data.at[index,"average_relative_ASA"] = round(statistics.mean([dssp[(model_chain.get_id(), residue_id)][3] for residue_id in matched_residue_ids]),2)
        threshold = 0.25
        data.at[index,"relative_ASA_below_025"] = round(sum([1 for residue_id in matched_residue_ids if dssp[(model_chain.get_id(), residue_id)][3] <= threshold])/len(matched_residue_ids),3)
        


  
    # save output 
    data.to_csv(csv_out, index = False, sep=',', header=True, na_rep='NA')



    
if __name__ == '__main__':
    main()