#!/bin/bash


for i in {1..22}
do
	for p in AFA HIS CHN
	do
		echo $i $p
		tail -n +42 /home/ryan/software/ensembl-vep/${p}_anno/protein_altering_snps_chr${i}.txt > /home/ryan/TOPMed_Proteome/06Elastic_net/PAV_modeling_methods/01PAV_annotation/PAV_annotations/${p}_protein_altering_snps_chr${i}.txt
	done
done
