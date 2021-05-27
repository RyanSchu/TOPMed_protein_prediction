###manhattan plot custom code adapted from https://www.r-graph-gallery.com/101_Manhattan_plot.html

library(ggplot2)
library(dplyr)
library(data.table)
library(Rmpfr)
library(argparse)
library(viridis)
library(qvalue)
library(tidyr)

parser <- ArgumentParser()
parser$add_argument("--cis", help="cis output from determine cis trans script")
parser$add_argument("--trans", help="trans output from determine cis trans script")
parser$add_argument("--title", help="plot title")
parser$add_argument("--out", help="prefix as in plink")
args <- parser$parse_args()

"%&%" = function(a,b) paste(a,b,sep="")

cis<-fread("zcat " %&% args$cis) %>% 
  separate(gene,into=c("aptamer","gene"),sep="_") 
trans<-fread("zcat " %&%  args$trans)

cis$cis_trans<-"cis"
trans$cis_trans<-"trans"
genome<-rbind.data.frame(cis,trans)
f<-qvalue(genome$pvalue)
genome<-genome %>% 
  filter(pvalue < 1e-4)
genome<- genome %>% separate(snps,into=c("chr","pos"),remove = F) %>% mutate(chr=as.numeric(gsub("chr","",chr)),pos=as.numeric(pos),cis_trans=ifelse(FDR>0.05,"FDR>0.05",cis_trans))
str(genome)
final_genome<-genome %>% group_by(chr) %>%
  summarise(chr_len=max(pos)) %>% 
  mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>% 
  select(-chr_len) %>%
  left_join(genome,.,by="chr") %>%
  arrange(chr,pos) %>%
  mutate(poscum=pos+tot)
axisdf<- final_genome %>% group_by(chr) %>% summarize(center=(sum(c(max(poscum) ,min(poscum))))/2)

if (args$title == "CAU"){
  args$title <- "EUR"
}

png(args$out %&% ".png",width = 1600,height=800,)
  ggplot(final_genome,aes(poscum,y=-log10(pvalue))) +
    geom_point(aes(color=as.factor(cis_trans)),alpha=0.8,size=1.3)+
    scale_x_continuous(label = axisdf$chr,breaks = axisdf$center) +
    theme_bw(25) +
    theme( 
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )+
    xlab("Chr") +
    ggtitle(args$title %&% " Cis/Trans Proteome QTLs") + 
    scale_color_manual(values=c("purple","black","orange"))
dev.off()
