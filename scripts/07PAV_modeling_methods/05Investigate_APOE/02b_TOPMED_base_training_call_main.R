setwd("/home/rschubert1/scratch/TOPMed_Proteome/06Elastic_net/PAV_modeling")
source("01a_EN_baseline.R")
"%&%" <- function(a,b) paste(a,b, sep='')

argv <- commandArgs(trailingOnly = TRUE)
pop <- argv[1]
chrom <- 19

expr="/data/rschubert1/TOPMED_Proteome/" %&% pop %&% "/PCAIR_modeling/03adjusted_expression/Proteome_TOPMed_" %&% pop %&% "_ln_adjAgeSex_mean_rank-inverse_adj10PCair_MEQTL_reformatted_ids.txt.gz"


snp_annot_file <- "/scratch/rschubert1/TOPMed_Proteome/investigate_APOE/APOE_genotypes/predd_input/uniq_pred_db_" %&% pop %&% "_hg38.chr19.maf0.01.R20.8.anno.txt.gz"
gene_annot_file <- "/scratch/rschubert1/TOPMed_Proteome/investigate_APOE/model_PAV_adj_APOE/APOE_annotation.txt"
genotype_file <- "/scratch/rschubert1/TOPMed_Proteome/investigate_APOE/APOE_genotypes/predd_input/uniq_pred_db_" %&% pop %&% "_hg38.chr19.maf0.01.R20.8.geno.txt.gz"
expression_file <- expr
prefix <- "/home/rschubert1/scratch/TOPMed_Proteome/investigate_APOE/model_PAV_adj_APOE/out/" %&% pop %&% "_chr" %&% chrom %&% "_base_APOE"
pav_map_file<-"/scratch/rschubert1/TOPMed_Proteome/06Elastic_net/PAV_modeling_methods/01PAV_annotation/PAV_annotations/ALL_protein_altering_snps_chr" %&% chrom %&% ".txt"
aptamer_ens_map<-"/scratch/rschubert1/TOPMed_Proteome/06Elastic_net/PAV_modeling_methods/01PAV_annotation/Soma_to_Transcript_mapping.csv"

main(snp_annot_file, gene_annot_file, genotype_file, expression_file, covariates_file, as.numeric(chrom), prefix, null_testing=FALSE,PAV_mapping_file=pav_map_file,aptamer_ensemble_annotation_file=aptamer_ens_map)


