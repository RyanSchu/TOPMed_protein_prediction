library(GENESIS)
library(SNPRelate)
library(GWASTools)
library(dplyr)
library(tibble)
library(data.table)
##LD prune

args = commandArgs(trailingOnly=TRUE)
"%&%" = function(a,b) paste(a,b,sep="")

gdsfile<-"/home/ryan/topmed/proteome/" %&% args[1] %&%  "/PCAIR/00autosome.gds"
# gds <- snpgdsOpen(gdsfile)
# snpset <- snpgdsLDpruning(gds, method="corr", slide.max.bp=10e6,
#                           ld.threshold=sqrt(0.3), verbose=FALSE)
# pruned <- unlist(snpset, use.names=FALSE)

# saveRDS(pruned,"/home/ryan/topmed/proteome/ALL/unrelated_ALL/PCair_steps/pruned_set.RDS")

# #calculate the kinship coefficients
# king<-snpgdsIBDKING(gds)
# kingMat<-king$kinship
# snpgdsClose(gds)
# 
# saveRDS(king,file="/home/ryan/topmed/proteome/ALL/unrelated_ALL/PCair_steps/King_matrix.RDS")
pruned<-readRDS("/home/ryan/topmed/proteome/" %&% args[1] %&%  "/PCAIR/pruned_set.RDS")
king<-readRDS("/home/ryan/topmed/proteome/" %&% args[1] %&%  "/PCAIR/King_matrix.RDS")
kingMat<-king$kinship
colnames(kingMat)<-king$sample.id
row.names(kingMat)<-king$sample.id
#create genotype object
geno <- GdsGenotypeReader(filename = gdsfile)
genoData <- GenotypeData(geno)

#run PC air

mypcair <- pcair(genoData, kinobj = kingMat, divobj = kingMat,
                 snp.include = pruned)
png("/home/ryan/topmed/proteome/" %&% args[1] %&%"/PCAIR/PCAIR_PC1_VS_PC2.png")
plot(mypcair)
dev.off()

eigenvec<-mypcair$vectors %>% as.data.frame() %>% rownames_to_column(var="sample_id")
str(eigenvec)
val<-mypcair$values %>% as.data.frame()
fwrite(eigenvec,"/home/ryan/topmed/proteome/" %&% args[1] %&%  "/PCAIR/PCAIR.eigenvec",col.names = T,row.names = F,sep='\t')
fwrite(val,"/home/ryan/topmed/proteome/" %&% args[1] %&%  "/PCAIR/PCAIR.eigenval")
