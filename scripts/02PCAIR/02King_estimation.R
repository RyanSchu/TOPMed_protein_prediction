library(GENESIS)
library(SNPRelate)
library(GWASTools)
library(dplyr)
library(tibble)

args = commandArgs(trailingOnly=TRUE)
"%&%" = function(a,b) paste(a,b,sep="")

##LD prune
gdsfile<-"/home/ryan/topmed/proteome/" %&% args[1] %&%  "/PCAIR/00autosome.gds"
gds <- snpgdsOpen(gdsfile)
# snpset <- snpgdsLDpruning(gds, method="corr", slide.max.bp=10e6, 
#                           ld.threshold=sqrt(0.3), verbose=FALSE)
# pruned <- unlist(snpset, use.names=FALSE)

#calculate the kinship coefficients
king<-snpgdsIBDKING(gds)
kingMat<-king$kinship
snpgdsClose(gds)

saveRDS(king,file="/home/ryan/topmed/proteome/" %&% args[1] %&%  "/PCAIR/King_matrix.RDS")
# #create genotype object
# geno <- GdsGenotypeReader(filename = gdsfile)
# genoData <- GenotypeData(geno)
# 
# #run PC air
# 
# mypcair <- pcair(genoData, kinobj = kingMat, divobj = kingMat,
#                  snp.include = pruned)
# png("/home/ryan/topmed/proteome/ALL/unrelated_ALL/PCair_steps/PCAIR_PC1_VS_PC2.png")
# plot(mypcair)
# dev.off()
# 
# eigenvec<-mypcair$vectors %>% rownames_to_column(var="sample_id")
# 
# fwrite(eigenvec,"/home/ryan/topmed/proteome/ALL/unrelated_ALL/PCair_steps/PCAIR.eigenvec",col.names = T,row.names = F,sep='\t')
# fwrite(mypcair$values,"/home/ryan/topmed/proteome/ALL/unrelated_ALL/PCair_steps/PCAIR.eigenval")
