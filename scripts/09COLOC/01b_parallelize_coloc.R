library(dplyr)
library(data.table)
source('/home/rschubert1/scratch/test_coloc_for_elyse/05_coloc.R')
"%&%" = function(a,b) paste(a,b,sep="")
# should have an input list of the genes to test in which populations and which chrom they are on?

coloc_analysis <- function(gene_id=NULL, pop=NULL, pop_size=NULL, phenotype=NULL){
	gene <- gsub("\\..*","",gene_id)
	print(gene)
	rds<-'/home/rschubert1/scratch/test_coloc_for_elyse/coloc/gene_list_' %&% pop %&% "_" %&% phenotype %&% ".RDS"
  saveRDS(gene_id,rds )
	eqtl <- paste('/home/egeoffroy/LD_matrix/coloc/pQTL_', pop, '_', phenotype, '.txt.gz', sep = '')
	gwas <- paste('/home/egeoffroy/LD_matrix/coloc/GWAS_TOPMED', pop, '_', phenotype, '.txt.gz', sep = '')
	
	if(!file.exists(eqtl)){
	  cat(eqtl, "not exists. Exiting.\n")
	  return(NA)
	}
	if(!file.exists(gwas)){
	  cat(eqtl, "not exists. Exiting.\n")
	  return(NA)
	}
	
        ld_matrix <- 'T'
		if(!is.null(ld_matrix)){
					F_gwas <- fread(gwas, header = T, stringsAsFactors=F)
					# print(F_gwas)
					gwas_size <- F_gwas$`sample_size`[1]
					# print(gwas_size)
					if(pop == 'CAU'){ # for some reason CAU uses pvalues instead of bse values
#						main(eqtl=eqtl, gwas=gwas, directory=paste('/home/egeoffroy/LD_matrix/', pop, '_1Mb_coords_LDMatrix', sep =''),  mode = 'p', gene_list=gene_id, eqtlGeneCol='gene_id', eqtlSNPCol='variant_id', eqtlMAFCol='maf', eqtlSeCol=6, eqtlBetaCol=5, eqtlSampleSize= pop_size, gwasSNPCol=1, LD=ld_matrix, method="cond", gwasBetaCol=2, gwasSeCol=3, gwasSampleSize=gwas_size, outFile=paste('/home/egeoffroy/LD_matrix/coloc_output/', pop, '/', pop, '_', phenotype, '.txt', sep = ''), ld_header = 'T')
						# eqtl1 <- fread(eqtl, header = T, stringsAsFactors=F)
						# gwas1 <- fread(gwas, header = T, stringsAsFactors=F)
						# eqtl1$variant_id <- str_replace_all(eqtl1$variant_id, 'chr', '')
						# gwas1$panel_variant_id <- str_replace_all(gwas1$panel_variant_id, 'chr', '')
						# write.table(eqtl1, eqtl, quote= F, row.names=F)
						# write.table(gwas1, gwas, quote=F, row.names=F)
					  
					  # eqtl <- paste('/home/rschubert1/scratch/test_coloc_for_elyse/coloc/pQTL_', pop, '_', phenotype, '.txt.gz', sep = '')
					  # gwas <- paste('/home/rschubert1/scratch/test_coloc_for_elyse/coloc/GWAS_TOPMED', pop, '_', phenotype, '.txt.gz', sep = '')
						LD_dir<-"/home/rschubert1/scratch/test_coloc_for_elyse/test_LD/CAU_1Mb_coords_LDMatrix"
					}  else {
					  LD_dir<-paste('/home/egeoffroy/LD_matrix/', pop, '_1Mb_coords_LDMatrix', sep ='')
					}
				  	command<-"qsub -v pop=" %&% pop %&% 
				  	         ",phenotype=" %&% phenotype %&%
				  	         ",eqtl=" %&% eqtl %&%
				  	         ",gwas=" %&% gwas %&%
				  	         ",LD_dir=" %&% LD_dir %&%
				  	         ",gene_id=" %&% rds %&%
				  	         ",pop_size=" %&% pop_size %&%
				  	         ",ld_matrix=" %&% ld_matrix %&%
				  	         ",gwas_size=" %&% gwas_size %&%
				  	         " -N " %&% pop %&% "_" %&% phenotype %&% "_coloc " %&%
				  	         "/scratch/rschubert1/test_coloc_for_elyse/qsub_wrapper.sh"
				  	system(command)
				
				}
}
