# GPCR_project

This directory contains the data and scripts used for the structural analysis in the context of *publication*.

## METHODS AND MATERIALS

1) The raw input data is stored under *data/221107_raw_data.csv*. We replace the entry 
```
253,delta,HCAR2;HCAR3,Q8TDS4;P49019,ENSG00000182782;ENSG00000255398,TRUE,TRUE,FALSE,,51,NRCLQRKMTGEPDNNRSTSVELTGDPNKTRGAPEALMANSGEPWSPSYLGP
```
with 
```
253,delta,HCAR2,Q8TDS4,ENSG00000182782,TRUE,TRUE,FALSE,,51,NRCLQRKMTGEPDNNRSTSVELTGDPNKTRGAPEALMANSGEPWSPSYLGP

253,delta,HCAR3,P49019,ENSG00000255398,TRUE,TRUE,FALSE,,51,NRCLQRKMTGEPDNNRSTSVELTGDPNKTRGAPEALMANSGEPWSPSYLGP
```
to ensure that we have a unique uniprot id in every row. We store the resulting data as *data.csv* which will be used in the subsequent analysis. The processed data is stored under *analysis/data_results.csv*.

2) From *https://alphafold.ebi.ac.uk/download* we download all structures under *human proteom* (*--2022-11-11 15:08:40--  https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/UP000005640_9606_HUMAN_v4.tar*). We then select the relevant models by uniprot id and store them under data/model_pdb.
3) Next we compute an alignment of the antigen sequence (*Antigen_Sequence*) with the protein sequence from the uniprot id (*Uniprot_ID*). When computing the alignment, identical characters are given 1 points, 0.5 points are deducted for each non-identical character, 1 points are deducted when opening a gap, and 0.5 points are deducted when extending it. We normalize this score by the length of the antigen sequence and save it as *normalized_alignment_score*. In particular, in case of a perfect match the score is 1. We also save the aligned subsequence of the protein as *aligned_subsequence*.
4) Next we analyze the region in the alphafold pdb file defined by the *aligned_subsequence* computed in the previous step. We compute:
    + the average plDDT score (*average_plDDT_score*) (that is, we sum up the plDDT scores for all residues in the subsequence and divide the result by the number of residues in the subsequence)
    + the number of plDDT scores below 50.0 divided by the total number of residues in the aligned subsequence (*plDDT_below_50*)
    + the average relative ASA (accessible solvent area) score (*average_relative_ASA*) (that is, we sum up the relative ASA scores for all residues in the subsequence and divide the result by the number of residues in the subsequence)
    + the number of relative ASA scores below 0.25 divided by the total number of residues in the aligned subsequence (*relative_ASA_below_025*)
    + the secondary structure (by 3-state characters, see DSSP below) (*secondary_structure*)
    + the number of states H in the secodary structure divided by the the total number of residues in the aligned subsequence (*relative_H-dssp*) (and similar scores for the states E and C)
5) For the final analysis and for producing figures we only use entries with a normalized_alignment_score greater or equal to 0.90 (note that a score of 1 means a perfect match)


# DSSP
+ we choose accessible surface area (ASA) values according to Sander and Rost, 1994, https://doi.org/10.1002/prot.340200303
+ see also https://biopython.org/docs/1.75/api/Bio.PDB.DSSP.html
+ for installation of dssp see: https://ssbio.readthedocs.io/en/latest/instructions/dssp.html 

| 8-letter character | 3-state character | Structure                    |
|--------------------|-------------------|------------------------------|
| H                  | H                 | Alpha helix (3-12)           |
| B                  | E                 | Isolated beta-bridge residue |
| E                  | E                 | Strand                       |
| G                  | H                 | 3-10 helix                   |
| I                  | H                 | Pi helix                     |
| T                  | C                 | Turn                         |
| S                  | C                 | Bend                         |
| -                  | C                 | None                         |


[comment]: # '## CMDs'

[comment]: # 'create list of all uniprot_ids: ```cat data.csv | cut -d ',' -f4 | sed '1d' | sed 's/;/\n/g' | sort | uniq > uniprot_ids.txt```'

