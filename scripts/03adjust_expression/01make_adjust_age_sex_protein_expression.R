## 00 make ALL data adjusted for age sex and exam

#################################################
# SET UP ENVIRONMENT
#################################################

library(dplyr)
library(ggplot2)
library(readxl)
library(data.table)
library(viridis)

"%&%" = function(a,b) paste(a,b,sep="")

data.dir = "/home/wheelerlab3/Data/TOPMed/"
gt.dir = "/home/wheelerlab3/topmed/gt_QC/"
prot.dir = data.dir %&% "TOPMed_Proteomics_MESA/"
out.dir = "/home/ryan/TOPMed_Proteome/03adjust_expression/"

################################################
# READ SAMPLE INFO 
################################################

info <- read_excel(data.dir %&% "TOPMed_ID_mapping/MESA_TOPMed_WideID_20190517.xlsx",skip = 2)
afa <- fread(gt.dir %&% "AFA/chr22.maf01.r2_8.psam")
chn <- fread(gt.dir %&% "CHN/chr22.maf01.r2_8.psam")
his <- fread(gt.dir %&% "HIS/chr22.maf01.r2_8.psam")
#get CAU info from Ryan
ryan.info <- fread("/home/ryan/topmed/proteome/sample_lists/proteome_basic_pop_codes.txt")
cau <- dplyr::filter(ryan.info, pop=="CAU") %>% select(sidno, gender_code) %>% rename(`#IID`=sidno, SEX=gender_code)
cau <- cau[duplicated(cau)==FALSE,]

################################################
# READ PROTEIN DATA
################################################

prot <- read_excel(prot.dir %&% "proteomics_data_merged_with_runlist_key_updated_mar1.xlsx", skip=7)
#FLAG in column RowCheck should be removed, see Proteomics_samples_to_exclude.xlsx
#rm FLAG rows
prot <- dplyr::filter(prot, RowCheck != "FLAG")


#pull relavent columns from info
#top_idX are the proteome sample ID's and top_tom_idX are the metabolome ID's for the same sample
protinfo <- dplyr::select(info,sidno,age1c,age5c,top_tom_id1,top_id1,top_ageatdraw1,top_tom_id5,top_id5,top_ageatdraw5)

#add exam variable to metab specifying exam 1 or exam 5 and add age
newprot <- mutate(prot, exam=ifelse(TOP_ID %in% protinfo$top_id1,1,NA))
prot <- mutate(newprot, exam=ifelse(TOP_ID %in% protinfo$top_id5,5,exam))

exam1prot <- filter(prot, exam==1)
exam5prot <- filter(prot, exam==5)

################################################
# JOIN WITH DEMOGRAPHIC INFO
################################################

protdf <- inner_join(prot, ryan.info, by=c('sidno','exam','TOP_ID'))
protdf <- dplyr::select(protdf, sidno, exam, age, gender_code, pop, starts_with("SL"), starts_with("HCE")) %>%
  select(-SlideId)


afa1df <- filter(protdf,pop=="AFA",exam==1,complete.cases(protdf))
afa5df <- filter(protdf,pop=="AFA",exam==5,complete.cases(protdf))
cau1df <- filter(protdf,pop=="CAU",exam==1,complete.cases(protdf))
cau5df <- filter(protdf,pop=="CAU",exam==5,complete.cases(protdf))
chn1df <- filter(protdf,pop=="CHN",exam==1,complete.cases(protdf))
chn5df <- filter(protdf,pop=="CHN",exam==5,complete.cases(protdf))
his1df <- filter(protdf,pop=="HIS",exam==1,complete.cases(protdf))
his5df <- filter(protdf,pop=="HIS",exam==5,complete.cases(protdf))
all1df <- filter(protdf,exam==1,complete.cases(protdf)); #print(all1df[!complete.cases(all1df),1])
all5df <- filter(protdf,exam==5,complete.cases(protdf)); #print(all5df[!complete.cases(all5df),1])


###############################################
# ADJUST
###############################################

library(preprocessCore) #has normalize.quantiles function
rankinv <- function(x){qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))}

adjmatlist <- list() #list to store adjmat's
adjdflist <- list() #list to store adjdf's
for(exam in c("1","5")){
  for ( pop in c("AFA","CHN","HIS","ALL","CAU"))
  {
    df <- get(tolower(pop) %&% exam %&% "df")
    rawmat <- as.matrix(df[,6:ncol(df)]); #print(dim(rawmat))
    logmat <- log(rawmat); #print(dim(logmat)) #natural log transform
    lmfunc <- function(x){resid(lm(x ~ df$age + df$gender_code))} #get residuals of prot after adj for age & sex
    adjmat <- apply(logmat, 2, lmfunc) ; #print(dim(adjmat)) #apply lmfunc to each column of logmat
    name <- pop %&% exam
    adjmatlist[[name]] <- adjmat
    adjdf <- cbind(df[,1], adjmat)
    adjdflist[[name]] <- adjdf
  }
}

elementMean <- function(my.list) { 
  arr <- array(unlist(my.list),c(dim(my.list[[1]])[1],dim(my.list[[1]])[2],length(my.list)))
  rowMeans( arr , na.rm=TRUE, dims = 2 )
}

# print(str(adjdflist))
#full join to add NA's to df if missing Exam 1 or 5
afadf <- full_join(adjdflist$AFA1, adjdflist$AFA5, by = "sidno"); print("afa")
caudf <- full_join(adjdflist$CAU1, adjdflist$CAU5, by = "sidno"); print("cau")
chndf <- full_join(adjdflist$CHN1, adjdflist$CHN5, by = "sidno"); print("chn")
hisdf <- full_join(adjdflist$HIS1, adjdflist$HIS5, by = "sidno"); print("his")
alldf <- full_join(adjdflist$ALL1, adjdflist$ALL5, by = "sidno"); print("all")

#################################
# FINAL ADJUSTMENTS AND WRITE
#################################

for(pop in c('AFA','CAU','CHN','HIS','ALL')){
  df <- get(tolower(pop) %&% "df")
  df1na <- select(df,ends_with(".x"))
  df5na <- select(df,ends_with(".y"))
  meanmat <- elementMean(list(df1na, df5na)) ##take the mean of exam 1 and exam 5
  invmeanmat <- round(apply(meanmat, 2, rankinv),6) ##rank-inverse normalize and round
  finaldf <- cbind(df[,1],as.data.frame(invmeanmat)) ##add sidno to df
  colnames(finaldf) <- c("sidno", colnames(protdf[6:ncol(protdf)])) ##retrieve column names
  fwrite(finaldf, file = out.dir %&% "Proteome_TOPMed_" %&% pop %&% "_ln_adjAgeSex_mean_rank-inverse.txt",quote=F,sep="\t")
}
