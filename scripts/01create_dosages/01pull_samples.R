library(dplyr)
library(readxl)
library(tidyr)
library(data.table)

#Read in relevant topmed ids
topmed<-read_xlsx("Z:/topmed/MESA_TOPMed_WideID_20190517.xlsx",skip = 2) 
names(topmed)
topmed_ids<-dplyr::select(topmed,sidno,top_id1,top_id5,top_ageatdraw1,top_ageatdraw5)

#separate by exam, process each, and determine n at both timepoints
meta_exam1<-select(topmed_ids,sidno,top_id1,top_ageatdraw1)
meta_exam1<-meta_exam1[complete.cases(meta_exam1),]
meta_exam1$exam<-1
colnames(meta_exam1)<-c("sidno","top_id","age","exam")

meta_exam5<-select(topmed_ids,sidno,top_id5,top_ageatdraw5)
meta_exam5<-meta_exam5[complete.cases(meta_exam5),]
meta_exam5$exam<-5
colnames(meta_exam5)<-c("sidno","top_id","age","exam")
meta_all<-rbind.data.frame(meta_exam1,meta_exam5)

proteome_ids<-read_xlsx("Z:/topmed/xlxs/proteome/proteomics_data_merged_with_runlist_key_updated_mar1.xlsx",skip=7) %>% select(sidno,TOP_ID)
proteome_with_age<-left_join(proteome_ids,meta_all,by=c("sidno","TOP_ID"="top_id"))

AFA_topmed_genotyped<-fread("Z:/topmed/gwas_qc/AFA/22/missingness_hwe_steps/05filtered_HWE.fam",drop = c("V2","V3","V4","V5","V6"))
AFA_topmed_genotyped$pop<-"AFA"
CHN_topmed_genotyped<-fread("Z:/topmed/gwas_qc/CHN/22/missingness_hwe_steps/05filtered_HWE.fam",drop = c("V2","V3","V4","V5","V6"))
CHN_topmed_genotyped$pop<-"CHN"
HIS_topmed_genotyped<-fread("Z:/topmed/gwas_qc/HIS/22/missingness_hwe_steps/05filtered_HWE.fam",drop = c("V2","V3","V4","V5","V6"))
HIS_topmed_genotyped$pop<-"HIS"

topmed_genotyped<-rbind.data.frame(AFA_topmed_genotyped,CHN_topmed_genotyped) %>% rbind.data.frame(HIS_topmed_genotyped)
topmed_genotyped$source<-"topmed"
colnames(topmed_genotyped)<-c("sidno","pop","source")

prot_with_age_pop<-left_join(proteome_with_age,topmed_genotyped,by="sidno")
#prot_with_age_no_pop<-filter(prot_with_age_pop, is.na(pop))

old_mesa_genotyped<-fread("Z:/mesa_files/dbgap_files/v6_SHARe_merged_groups/combined_consent_basic_pop_codes.txt")
old_mesa_genotyped$source<-"mesa"

prot_topmed_mesa<-left_join(prot_with_age_pop,old_mesa_genotyped,by="sidno") %>% mutate(pop.x = ifelse(is.na(pop.x), pop.y,pop.x)) %>% mutate(source.x = ifelse(is.na(source.x), source.y,source.x))
names(prot_topmed_mesa)
prot_topmed_mesa<-select(prot_topmed_mesa,sidno,TOP_ID,age,exam,pop.x,source.x,dbGaP_Subject_ID,pop_code,gender_coide)
colnames(prot_topmed_mesa)<-c("sidno","TOP_ID","age","exam","pop","source","dbGaP_Subject_ID","pop_code","gender_code")

proteome_with_age<-proteome_with_age[!is.na(proteome_with_age$TOP_ID),]
cat("Proteome n in each pop Exam 1\n")
for(i in c('AFA','CHN','HIS','CAU')){
  #geno <- get(tolower(i))
  n = dim(filter(prot_topmed_mesa , pop==i,exam ==5))[1]
  cat(i, ": ", n, "\n", sep="")
}

cat("Proteome n in each pop Exam 5\n")
for(pop in c('AFA','CHN','HIS','CAU')){
  geno <- get(tolower(pop))
  n = dim(filter(exam5prot,sidno %in% geno$`#IID`))[1]
  cat(pop, ": ", n, "\n", sep="")
}