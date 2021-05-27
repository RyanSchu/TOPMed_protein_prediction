# Author: Elyse Geoffroy
# This script calls function in 05b to run coloc. Parameter required is population id 
# Note that population sample sizes are hard-coded for the MESA populations. This will have to become a new parameter in COMP BIO project pipeline
library(dplyr)
library(stringr)
library(data.table)
source('/home/rschubert1/scratch/test_coloc_for_elyse/05b_run_coloc.R')

# Get command line arguments
args <- commandArgs(trailingOnly=TRUE)
print(args)
pop <- args[1]
print(pop)


files <- list.files('/home/rschubert1/scratch/test_coloc_for_elyse/SPred_results/bonf_filtered_results/', pattern = '^PCAIR', recursive = F, full.names=T)
SPrediXcan_Results <- data.frame()
for(file in files){
	file <- fread(file, header= T, sep = ',', stringsAsFactors=F)
	SPrediXcan_Results <- rbind(SPrediXcan_Results, file)
}
#print(head(SPrediXcan_Results))

#pops <- c('AFA')
#pops <- c('AFA', 'ALL', 'CAU', 'CHN', 'HIS')

if(pop=='AFA'){sample_size<-183}
if(pop=='ALL'){sample_size<-971}
if(pop=='CAU'){sample_size<-416}
if(pop=='CHN'){sample_size<-71}
if(pop=='HIS'){sample_size<-301}

SPred_pop <- SPrediXcan_Results %>% filter(Model == pop) %>% select(gene, Phenotype)
SPred_pop$Phenotype <- str_replace_all(SPred_pop$Phenotype, 'LDL_choleseterol', 'LDL_cholesterol')
# print(SPred_pop)

#for(i in 1:nrow(SPred_pop)){
#	row <- SPred_pop[i,]
#	coloc_analysis(pop=as.character(pop), pop_size=sample_size, phenotype=row$Phenotype, gene_id=row$gene)  
#}
phenos <- unique(SPred_pop$Phenotype)
for(pheno in phenos){
	data <- SPred_pop %>% filter(Phenotype == pheno)
	genes <- unique(data$gene)
	coloc_analysis(pop=as.character(pop), pop_size=sample_size, phenotype=pheno, gene_id=genes)
  
}
