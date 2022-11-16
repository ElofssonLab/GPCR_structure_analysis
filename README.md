# GPCR_project

## METHODS - single structure with single protein as antigen

1) we replace the entry 
```
253,delta,HCAR2;HCAR3,Q8TDS4;P49019,ENSG00000182782;ENSG00000255398,TRUE,TRUE,FALSE,,51,NRCLQRKMTGEPDNNRSTSVELTGDPNKTRGAPEALMANSGEPWSPSYLGP
```
with 
```
253,delta,HCAR2,Q8TDS4,ENSG00000182782,TRUE,TRUE,FALSE,,51,NRCLQRKMTGEPDNNRSTSVELTGDPNKTRGAPEALMANSGEPWSPSYLGP

253,delta,HCAR3,P49019,ENSG00000255398,TRUE,TRUE,FALSE,,51,NRCLQRKMTGEPDNNRSTSVELTGDPNKTRGAPEALMANSGEPWSPSYLGP
```
to ensure that we have a unique uniprot id in every row.

2) From *https://alphafold.ebi.ac.uk/download* we download all structures under *human proteom* (*--2022-11-11 15:08:40--  https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/UP000005640_9606_HUMAN_v4.tar*)
3) Next we compute an alignment of the antigen sequence with the protein sequence from the uniprot id. In the alignment score, identical characters are given 1 points, 0.5 point is deducted for each non-identical character. 1 points are deducted when opening a gap, and 0.5 points are deducted when extending it. We normalize this score by the length of the antigen sequence and save it as *normalized_alignment_score*. In particular, in case of a perfect match the score is 1. We also save the aligned subsequence of the protein as *aligned_subsequence*.
4) Next we analyze the region in the alphafold pdb file defined by the *aligned_subsequence* computed in the previous step. We compute:
    + the average plDDT score (*average_plDDT_score*)
    + the number of plDDT scores below 50.0 divided by the total number of residues in the aligned subsequence (*plDDT_below_50*)
    + the average relative ASA (accessible solvent area) score (*average_relative_ASA*)
    + the number of relative ASA scores below 0.25 divided by the total number of residues in the aligned subsequence (*relative_ASA_below_025*)
    + the secondary structure (by 3-state characters, see DSSP below) (*secondary_structure*)
    + the number of states H in the secodary structure divided by the the total number of residues in the aligned subsequence (*relative_H-dssp*) (and similar scores for the states E and C)



# DSSP
+ see also https://biopython.org/docs/1.75/api/Bio.PDB.DSSP.html
+ we choose accessible surface area (ASA) values according to Wilke; Tien et al. 2013 https://doi.org/10.1371/journal.pone.0080635

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

## CMDs

create list of all uniprot_ids: ```cat data.csv | cut -d ',' -f4 | sed '1d' | sed 's/;/\n/g' | sort | uniq > uniprot_ids.txt```


# Notes

See below link for installation of dssp:
https://ssbio.readthedocs.io/en/latest/instructions/dssp.html 




