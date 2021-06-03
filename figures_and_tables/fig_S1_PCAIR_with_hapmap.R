#Make nicely coloured pcaplots for the PCAIR PCs
library(data.table)
library(ggplot2)
library(viridisLite)
library(dplyr)
#library(cowplot)
"%&%" = function(a,b) paste (a,b,sep="")


hapmap_IDS<-fread("/home/ryan/HAPMAP3_hg19/hg38/topmed_id_format/pop_HM3_hg19_forPCA.txt",header=F) %>% 
  select(V1,V3)
colnames(hapmap_IDS)<-c("pop","sample_id")

for (population in c("AFA")){
  eigenval<-fread("/home/ryan/topmed/proteome/" %&% population %&% "/PCAIR_modeling/02PCAIR/PCAIR_with_hapmap_uncoloured.eigenval") %>% unlist() %>% unname()
  eval <- eigenval[1:10]
  scree<-round(eval/sum(eval),3)
  scree<-cbind.data.frame(scree,1:10)
  colnames(scree)<-c("percent_var","PC")
  eigenvec<-fread("/home/ryan/topmed/proteome/" %&% population %&% "/PCAIR_modeling/02PCAIR/PCAIR_with_hapmap_uncoloured.eigenvec")
  
  # print(str(hapmap_IDS))
  # print(str(eigenvec))
  pcdf<-full_join(hapmap_IDS,eigenvec,by="sample_id") %>% mutate(pop=if_else(is.na(pop),population,pop)) %>% mutate(pop=if_else(pop == "CAU","EUR",pop))

  pcaplots<-"/home/ryan/topmed/proteome/" %&% population %&% "/PCAIR_modeling/02PCAIR/PCAIR_with_hapmap_labels.pdf"  
  pdf(pcaplots)
  g1<-ggplot(data=scree, aes(x=PC, y=percent_var)) + 
    geom_point() + 
    geom_line() + 
    theme_bw() +
    scale_x_continuous(breaks = 1:10) +
    ggtitle("Proportion of variance explained") +
    xlab("Principal Component") +
    ylab("Percent Variance Explained")
  
  g1_alt<-ggplot(data=scree, aes(x=PC, y=percent_var)) + 
    geom_point() + 
    geom_line() + 
    theme_bw() +
    scale_x_continuous(breaks = 1:10) +
    xlab("Principal Component") +
    ylab("Percent Variance Explained")
  

  
  #PCA Plot 1 (PC1 vs PC2)
  g2<-ggplot() + 
    geom_point(data=pcdf,aes(x=V1,y=V2,col=pop,shape=pop)) + 
    theme_bw() + 
    scale_colour_brewer(palette="Set1") + 
    ggtitle("PC1 vs PC2") +
    xlab("PC1") +
    ylab("PC2")
  g2_alt<-ggplot() + 
    geom_point(data=pcdf,aes(x=V1,y=V2,col=pop,shape=pop)) + 
    theme_bw() + 
    scale_colour_brewer(palette="Set1") +
    xlab("PC1") +
    ylab("PC2")
  
  #PCA Plot 2 (PC1 vs PC3)
  g3<-ggplot() + 
    geom_point(data=pcdf,aes(x=V1,y=V3,col=pop,shape=pop)) + 
    theme_bw() + 
    scale_colour_brewer(palette="Set1") + 
    ggtitle("PC1 vs PC3") +
    xlab("PC1") +
    ylab("PC3")
  
  #PCA Plot 3 (PC2 vs PC3)
  g4<-ggplot() + 
    geom_point(data=pcdf,aes(x=V2,y=V3,col=pop,shape=pop)) + 
    theme_bw() + 
    scale_colour_brewer(palette="Set1") + 
    ggtitle("PC2 vs PC3") +
    xlab("PC2") +
    ylab("PC3")
  if (population == "ALL"){
    sample_list<-fread("/home/ryan/topmed/proteome/sample_lists/proteome_basic_pop_codes.txt",header=T) %>% 
      select(pop,sidno) %>% rename(sample_id=sidno) %>% rbind.data.frame(hapmap_IDS) %>% mutate(pop=if_else(pop == "CAU","EUR",pop))
    # str(sampl)s
    eigenvec<-mutate(eigenvec,sample_id=gsub("[0-9]+_","",sample_id))
    pcdf<-inner_join(sample_list,eigenvec,by="sample_id") 
    
    # print(str(pcdf))
    
    g2<-ggplot() + 
      geom_point(data=pcdf,aes(x=V1,y=V2,col=pop,shape=pop)) + 
      scale_shape_manual(values=1:7) +
      theme_bw() + 
      scale_colour_brewer(palette="Set1") + 
      ggtitle("PC1 vs PC2") +
      xlab("PC1") +
      ylab("PC2")
    
    g2_alt<-ggplot() + 
      geom_point(data=pcdf,aes(x=V1,y=V2,col=pop,shape=pop)) + 
      scale_shape_manual(values=1:7) +
      theme_bw() + 
      scale_colour_brewer(palette="Set1")  +
      xlab("PC1") +
      ylab("PC2")
    
    #PCA Plot 2 (PC1 vs PC3)
    g3<-ggplot() + 
      geom_point(data=pcdf,aes(x=V1,y=V3,col=pop,shape=pop)) + 
      scale_shape_manual(values=1:7) +
      theme_bw() + 
      scale_colour_brewer(palette="Set1") + 
      ggtitle("PC1 vs PC3") +
      xlab("PC1") +
      ylab("PC3")
    
    #PCA Plot 3 (PC2 vs PC3)
    g4<-ggplot() + 
      geom_point(data=pcdf,aes(x=V2,y=V3,col=pop,shape=pop)) + 
      scale_shape_manual(values=1:7) +
      theme_bw() + 
      scale_colour_brewer(palette="Set1") + 
      ggtitle("PC2 vs PC3") +
      xlab("PC2") +
      ylab("PC3")
  }
  if (population == "CAU"){
    population<-"EUR"
  }
  # print(str(g1_alt))
  g1_alt + ggsave("/home/ryan/TOPMed_Proteome/02PCAIR/hapmap_merge/screeplot_" %&% population %&% ".pdf",dpi=500)
  g1_alt + ggsave("/home/ryan/TOPMed_Proteome/02PCAIR/hapmap_merge/screeplot_" %&% population %&% ".png",dpi=500)
  g2_alt + ggsave("/home/ryan/TOPMed_Proteome/02PCAIR/hapmap_merge/PC1_vs_PC2_" %&% population %&% ".pdf",dpi=500)
  g2_alt + ggsave("/home/ryan/TOPMed_Proteome/02PCAIR/hapmap_merge/PC1_vs_PC2_" %&% population %&% ".png",dpi=500)
#  plot_grid(g2_alt,g1_alt) + ggsave("/home/ryan/TOPMed_Proteome/02PCAIR/hapmap_merge/scree_and_biplot_" %&% population %&% ".pdf",dpi=500)
#  plot_grid(g2_alt,g1_alt) + ggsave("/home/ryan/TOPMed_Proteome/02PCAIR/hapmap_merge/scree_and_biplot_" %&% population %&% ".png",dpi=500)
  dev.off()
}
#print(g1_alt)
#g1_