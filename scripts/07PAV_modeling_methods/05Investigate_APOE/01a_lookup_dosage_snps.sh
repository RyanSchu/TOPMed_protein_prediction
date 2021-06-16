#!/bin/bash

for i in AFA CHN HIS CAU ALL
do
	echo $i
	zgrep -f /scratch/rschubert1/TOPMed_Proteome/investigate_APOE/APOE_snps.txt /data/rschubert1/TOPMED_Proteome/${i}/PCAIR_modeling/01genotypes/dosages/no_header_hg38chr19.maf0.01.R20.8.dosage.txt.gz > /scratch/rschubert1/TOPMed_Proteome/investigate_APOE/${i}_dosage_snps.txt
done
