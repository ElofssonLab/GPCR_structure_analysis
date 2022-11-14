# GPCR_project

## QUESTIONS

Some entries have experimental structures


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
2) 

## CMDs

create list of all uniprot_ids: ```cat data.csv | cut -d ',' -f4 | sed '1d' | sed 's/;/\n/g' | sort | uniq > uniprot_ids.txt```

## TODO list

+ Download experimental pdb structures (97 available)
+ calculate alignment, get start and end residue number to define a region





# Notes

See below link for installation of dssp:
https://ssbio.readthedocs.io/en/latest/instructions/dssp.html 




# DSSP
see https://biopython.org/docs/1.75/api/Bio.PDB.DSSP.html

TODO short needs to be verified

| Code | Short | Structure                    |
|------|-------|------------------------------|
| H    | H     | Alpha helix (3-12)           |
| B    | E     | Isolated beta-bridge residue |
| E    | E     | Strand                       |
| G    | H     | 3-10 helix                   |
| I    | H     | Pi helix                     |
| T    | L     | Turn                         |
| S    | L     | Bend                         |
| -    | L     | None                         |

