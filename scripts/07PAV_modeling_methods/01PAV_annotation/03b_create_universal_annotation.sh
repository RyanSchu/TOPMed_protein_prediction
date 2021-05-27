#!/bin/bash


for i in {1..22}
do
	for p in CHN HIS CAU AFA
	do
		echo $p $i
		cut -f 2,4-7 -d$'\t' /home/ryan/TOPMed_Proteome/06Elastic_net/PAV_modeling_methods/01PAV_annotation/PAV_annotations/${p}_protein_altering_snps_chr${i}.txt >> /home/ryan/TOPMed_Proteome/06Elastic_net/PAV_modeling_methods/01PAV_annotation/PAV_annotations/ALL_protein_altering_snps_chr${i}.txt
	done
	sort -u -o /home/ryan/TOPMed_Proteome/06Elastic_net/PAV_modeling_methods/01PAV_annotation/PAV_annotations/ALL_protein_altering_snps_chr${i}.txt /home/ryan/TOPMed_Proteome/06Elastic_net/PAV_modeling_methods/01PAV_annotation/PAV_annotations/ALL_protein_altering_snps_chr${i}.txt
done
