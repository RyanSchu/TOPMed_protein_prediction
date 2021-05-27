setwd("/home/rschubert1/scratch/TOPMed_Proteome/06Elastic_net")
source("02a_dapg_nested_cv_elnet.R")
"%&%" <- function(a,b) paste(a,b, sep='')

argv <- commandArgs(trailingOnly = TRUE)
pop <- argv[1]
chrom <- argv[2]
pip <- argv[3]
clst <- argv[4]

if ( clst == "T" ) {
suffix="T"
filt=TRUE 
} else {
suffix="F"
filt=FALSE
}

if ( pop == "ALL"){
        expr="/data/rschubert1/TOPMED_Proteome/ALL/PCAIR_modeling/03adjusted_expression/Proteome_TOPMed_ALL_ln_adjAgeSex_mean_rank-inverse_adj10PCair_MEQTL_reformatted_ids_sorted.txt.gz"
} else {
        expr="/data/rschubert1/TOPMED_Proteome/" %&% pop %&% "/PCAIR_modeling/03adjusted_expression/Proteome_TOPMed_" %&% pop %&% "_ln_adjAgeSex_mean_rank-inverse_adj10PCair_MEQTL_reformatted_ids.txt.gz"
}

snp_annot_file <- "/data/rschubert1/TOPMED_Proteome/" %&% pop %&% "/PCAIR_modeling/01genotypes/preddb_input/uniq_pred_db_hg38.chr" %&% chrom %&% ".maf0.01.R20.8.anno.txt.gz"
gene_annot_file <- "/data/rschubert1/TOPMED_Proteome/annotation_all_aptamers_ENSG.txt"
genotype_file <- "/data/rschubert1/TOPMED_Proteome/" %&% pop %&% "/PCAIR_modeling/01genotypes/preddb_input/uniq_pred_db_hg38.chr" %&% chrom %&% ".maf0.01.R20.8.geno.txt.gz"
expression_file <- expr
dapg_file <- "/data/rschubert1/TOPMED_Proteome/" %&% pop %&% "/PCAIR_modeling/05dapg/summary_dapg_out/summary_snps.txt"
prefix <- "/home/rschubert1/data/TOPMED_Proteome/" %&% pop %&% "/PCAIR_modeling/06Elastic_net/dapg_" %&% pip %&% "_" %&% suffix %&% "_out/" %&% pop %&% "_chr" %&% chrom %&% "_" %&% "PIP_" %&% pip %&% "_clus_filt_" %&% suffix

main(snp_annot_file, gene_annot_file, genotype_file, expression_file, dapg_file,covariates_file, as.numeric(chrom), prefix, null_testing=FALSE,filt_cluster=filt, pip_threshold=pip)


