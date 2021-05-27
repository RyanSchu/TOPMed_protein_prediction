#!/bin/bash
directory='/home/egeoffroy/Wojcik/Wojcik_build38/gzipped_versions'
for file in `ls $directory`; do
        file='/home/egeoffroy/Wojcik/Wojcik_build38/gzipped_versions/'${file}
        name=$(echo ${file}| cut -d '/' -f 7)
        name=$(echo ${name}| cut -d'.' -f 1)
        suffix='WojcikG_'
        name=${name#"$suffix"}
        echo ${name}
        nohup bash /home/rschubert1/scratch/TOPMed_Proteome/08PWAS/Wojcik_SPrediXcan_cvModels_Meta2.sh ${file} ${name} > nohup4_${name}.out &
done
