library(data.table)
library(dplyr)
library(argparse)
library(tidyr)

parser <- ArgumentParser()
parser$add_argument("--proteins", help="pqtl output by MEQTL")
parser$add_argument("--annotation", help="cis snps determined by ")
parser$add_argument("--out", help="file you would like to output as")
args <- parser$parse_args()

"%&%" = function(a,b) paste(a,b,sep="")

proteins<-fread("zcat " %&% args$proteins,header=T)
annotation<-fread("zcat " %&% args$annotation,header = F)

colnames(annotation)<-c("joint_id","c","start","end")

new<-separate(annotation, col = "joint_id",into=c("aptamer","gene"),sep="_",remove=F)
str(new)
final<-left_join(new,proteins,by=c("aptamer"="id")) %>% select(-c,-start,-end,-aptamer,-gene)

fwrite(final,args$out,col.names = T,sep='\t')



