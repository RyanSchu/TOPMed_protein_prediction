library(GENESIS)
library(SNPRelate)
library(GWASTools)
library(dplyr)
library(tibble)
library(data.table)
##LD prune

args = commandArgs(trailingOnly=TRUE)
"%&%" = function(a,b) paste(a,b,sep="")

if (args[1] == "CAU"){
  gdsfile<-"/home/ryan/topmed/proteome/CAU/genotypes/QC_lauren/hapmap/07final_LDpruned.gds"
  fam<-"/home/ryan/HAPMAP3_hg19/topmed_id_format/HM3_ASN_CEU_YRI_Unrelated_hg19_noAmbig.fam"
} else {
  gdsfile<-"/home/ryan/topmed/proteome/" %&% args[1] %&% "/genotypes/QC/hapmap/07final_LDpruned.gds"
  fam<-"/home/ryan/HAPMAP3_hg19/hg38/topmed_id_format/hapmap.hg38.topmed.id.format.fam"
}

hapmap.fam<-fread(fam,header=F)
hapmap_IDS<-unname(unlist(hapmap.fam$V2))
if (args[1] == "ALL"){
  king<-readRDS("/home/ryan/topmed/proteome/ALL/PCAIR_ALL/01PCair_steps/King_matrix.RDS")
} else {
  
  king<-readRDS("/home/ryan/topmed/proteome/" %&% args[1] %&% "/PCAIR_modeling/02PCAIR/King_matrix.RDS")
}
kingMat<-king$kinship
colnames(kingMat)<-king$sample.id
row.names(kingMat)<-king$sample.id

build_down<-matrix(-1,nrow=length(hapmap_IDS),ncol = ncol(kingMat))
colnames(build_down)<-colnames(kingMat)
row.names(build_down)<-hapmap_IDS

kingMat<-rbind(kingMat,build_down)

build_right<-matrix(-1,nrow=nrow(kingMat),ncol = length(hapmap_IDS))
colnames(build_right)<-hapmap_IDS
row.names(build_right)<-row.names(kingMat)

kingMat<-cbind(kingMat,build_right)


geno <- GdsGenotypeReader(filename = gdsfile)
genoData <- GenotypeData(geno)

#run PC air


mypcair <- pcair(genoData, kinobj = kingMat, divobj = kingMat,unrel.set=hapmap_IDS)
png("/home/ryan/topmed/proteome/" %&% args[1] %&% "/PCAIR_modeling/02PCAIR/PCAIR_with_hapmap_uncoloured_PC1_VS_PC2.png")
plot(mypcair)
dev.off()

eigenvec<-mypcair$vectors %>% as.data.frame() %>% rownames_to_column(var="sample_id")
str(eigenvec)
val<-mypcair$values %>% as.data.frame()
fwrite(eigenvec,"/home/ryan/topmed/proteome/" %&% args[1] %&% "/PCAIR_modeling/02PCAIR/PCAIR_with_hapmap_uncoloured.eigenvec",col.names = T,row.names = F,sep='\t')
fwrite(val,"/home/ryan/topmed/proteome/" %&% args[1] %&% "/PCAIR_modeling/02PCAIR/PCAIR_with_hapmap_uncoloured.eigenval")
