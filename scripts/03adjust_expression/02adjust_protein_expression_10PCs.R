## adjust protein expression for each population
library(data.table)
library(dplyr)
library(tidyr)
library(gtools)
library(tibble)
"%&%" = function(a,b) paste (a,b,sep="")

prot.dir = "/home/ryan/TOPMed_Proteome/03adjust_expression/expression/"
# pc.dir = "/home/ryan/topmed/proteome/CAU/PCAIR_modeling/02PCAIR/" 
annot.dir = "/home/wheelerlab3/topmed/annot_files/"

# 
# CAU <- fread(prot.dir %&% "Proteome_TOPMed_CAU_ln_adjAgeSex_mean_rank-inverse.txt", header=TRUE)
# AFA <- fread(prot.dir %&% "Proteome_TOPMed_AFA_ln_adjAgeSex_mean_rank-inverse.txt", header=TRUE)
# HIS <- fread(prot.dir %&% "Proteome_TOPMed_HIS_ln_adjAgeSex_mean_rank-inverse.txt", header=TRUE)
# CHN <- fread(prot.dir %&% "Proteome_TOPMed_CHN_ln_adjAgeSex_mean_rank-inverse.txt", header=TRUE)
# ALL<-rbind.data.frame(CAU,CHN) %>% rbind.data.frame(HIS) %>% rbind.data.frame(AFA)

for (pop in c("AFA","CHN","CAU","ALL","HIS")){
  print(pop)
  out.dir = "/home/ryan/topmed/proteome/" %&% pop %&% "/PCAIR_modeling/03adjusted_expression/"
  
  prot<-fread("zcat " %&% prot.dir %&% "Proteome_TOPMed_" %&% pop %&% "_ln_adjAgeSex_mean_rank-inverse.txt.gz", header=TRUE)
  # fwrite(prot, out.dir %&% "Proteome_TOPMed_" %&% pop %&% "_ln_adjAgeSex_mean_rank-inverse_0PCs.txt",sep='\t',quote=F)
  pcs <- fread("/home/ryan/topmed/proteome/" %&% pop %&% "/PCAIR_modeling/02PCAIR/PCAIR.eigenvec",header=TRUE)
  pcs<-mutate(pcs,sample_id=gsub("[0-9]+_","",sample_id))
  #make protein df in same order as PCs
  protdf <- filter(prot, sidno %in% pcs$sample_id) %>% arrange(sidno)
  pcsdf <- filter(pcs, sample_id %in% prot$sidno) %>% arrange(sample_id)
  stopifnot(protdf$sidno == pcsdf$sample_id)
  protmat <- as.matrix(protdf[,-1])
  pcmat <- as.matrix(pcsdf[,2:11])
  lmfunc <- function(x){resid(lm(x ~ pcmat))} #get residuals of prot after adj for 10PCs
  adjmat <- apply(protmat, 2, lmfunc) #apply lmfunc to each column of logmat
  rankinv <- function(x){qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))}
  invmat <- round(apply(adjmat, 2, rankinv),6) ##rank-inverse normalize and round
  #output invmat in MatrixEQTL format
  tprot <- t(invmat)
  colnames(tprot) <- as.character(protdf$sidno)
  a <- data.frame(tprot)
  b <- tibble::rownames_to_column(data.frame(a), "id")
  colnames(b) <- c("id",colnames(tprot))
  
  annot = fread(annot.dir %&% "somascan1.3k_gencode.v32_annotation.txt")
  
  annotinfo <- inner_join(b,annot,by=c("id"="SomaId"))
  #choose first position in annot as aptamer position
  annotinfo <- annotinfo[!duplicated(annotinfo$id),]
  #sort by chr/pos
  allsorted1 <- arrange(annotinfo, c, start) 
  allsorted <- allsorted1[gtools::mixedorder(allsorted1$c),] #sort chromosomes numerically, i.e. 1, 2, 3...X, Y not 1, 10, 11...X, Y
  proteome <- allsorted[1:ncol(b)]
  protloc <- select(allsorted,id,chr,start,end)
  colnames(protloc) <- c("id","chr","s1","s2")
  
  #27 aptamers didn't map to genome, write with NAs
  #I'm not sure if MatrixEQTL will allow it--we may have to run them separately
  fwrite(proteome, out.dir %&% "Proteome_TOPMed_" %&% pop %&% "_ln_adjAgeSex_mean_rank-inverse_adj10PCs_MEQTL.txt",sep='\t',na="NA",quote=F)
  fwrite(protloc, out.dir %&% "ProtLoc_TOPMed_" %&% pop %&% "_ln_adjAgeSex_mean_rank-inverse_adj10PCs_MEQTL.txt",sep='\t',na="NA",quote=F)
  
}




