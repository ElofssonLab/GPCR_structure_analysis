#!/bin/bash
set -e



git_directory=/home/samfro/git/GPCR_project

raw_files=$git_directory/data/human_proteome
data_csv=$git_directory/data/data.csv
model_pdb=$git_directory/data/model_pdb
model_cif=$git_directory/data/model_cif

# create ID
id_file=$git_directory/data/uniprot_ids.out
cat $data_csv | cut -d ',' -f4 | sed '1d' | sort | uniq > $id_file
# show file | select fourth column | delete first line |  sort | delete copies 

#cat $data_csv | cut -d ',' -f4 | sed '1d' | sed 's/;/\n/g' | sort | uniq > $id_file
## show file | select fourth column | delete first line | replace ";" with "\n" | sort | delete copies 



#id_file=$git_directory/data/uniprot_ids.txt
nber_lines=$(wc -l $id_file | awk '{print $1}')
echo $nber_lines


for n in $(seq 1 $nber_lines)
do 
        LN=$n
        id=$(sed -n ${LN}p $id_file) 

        pdb_path=$(find $raw_files -name "*${id}*\.pdb\.gz")
        cif_path=$(find $raw_files -name "*${id}*\.cif\.gz")

        echo $pdb_path
        echo $cif_path
        cp $pdb_path $model_pdb
        cp $cif_path $model_cif

done

gzip -d ${model_pdb}/*
gzip -d ${model_cif}/*