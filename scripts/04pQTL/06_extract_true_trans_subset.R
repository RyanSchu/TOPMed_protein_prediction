
library(data.table)
library(dplyr)
library(tidyr)

pop_list<-c("ALL")


"%&%" = function(a,b) paste(a,b,sep="")

annot<-fread("/home/ryan/topmed/proteome/somascan1.3k_gencode.v32_annotation.txt",header=T)

for (pop in pop_list){
  potential_trans<-fread("zcat /home/ryan/topmed/proteome/" %&% pop %&% "/PCAIR_modeling/04pQTL/04b_trans_QTLs/trans_eQTLs_" %&% pop %&% "_chr1_trans_1e-4.txt.gz") %>% 
    separate(gene,into=c("aptamer","gene"),sep="_")

  for (chr in 2:22) {
    cat(pop,chr,"\n")
    tmp<-fread("zcat /home/ryan/topmed/proteome/" %&% pop %&% "/PCAIR_modeling/04pQTL/04b_trans_QTLs/trans_eQTLs_" %&% pop %&% "_chr" %&% chr %&% "_trans_1e-4.txt.gz") %>% 
      separate(gene,into=c("aptamer","gene"),sep="_")
    potential_trans<-rbind.data.frame(potential_trans,tmp)
  }
  
  true_cis<-fread("zcat /home/ryan/topmed/proteome/" %&% pop %&% "/PCAIR_modeling/04pQTL/cis_eQTLs_" %&% pop %&% "_WG_all_cis.txt.gz") %>% 
    separate(gene,into=c("aptamer","gene"),sep="_")
  

  
  true_trans<-anti_join(potential_trans,true_cis,by=c("snps","aptamer"))

  fwrite(true_trans,"/home/ryan/topmed/proteome/" %&% pop %&% "/PCAIR_modeling/04pQTL/04b_trans_QTLs/" %&% pop %&% "true_trans_pQTL_1e-4.txt",sep='\t',col.names=T)
}
