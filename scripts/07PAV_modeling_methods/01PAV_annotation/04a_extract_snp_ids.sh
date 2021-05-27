#!/bin/bash

for i in {1..22}
do
	cut -f 1,1 -d$'\t' /home/ryan/TOPMed_Proteome/06Elastic_net/PAV_modeling_methods/01PAV_annotation/PAV_annotations/ALL_protein_altering_snps_chr${i}.txt >> PAV_list.txt
done
sort -u -o PAV_list.txt PAV_list.txt
