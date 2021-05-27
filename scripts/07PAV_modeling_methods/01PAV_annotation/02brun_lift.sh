#!/bin/bash

chr=$1

Rscript lift_cau_PAV_annotation.R \
/home/ryan/software/ensembl-vep/CAU_anno/protein_altering_snps_chr${chr}.txt \
/home/ryan/software/ensembl-vep/CAU_anno/liftover/lifted_hg38_chr${chr}.txt \
/home/ryan/TOPMed_Proteome/06Elastic_net/PAV_modeling_methods/01PAV_annotation/PAV_annotations/CAU_protein_altering_snps_chr${chr}.txt


