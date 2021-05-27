#lifotver cau
library(dplyr)
library(data.table)
"%&%" = function(a,b) paste0(a,b)

args <- commandArgs(trailingOnly = TRUE)
input<-args[1]
liftover<-args[2]
output<-args[3]

lift_map<-fread(liftover) %>%
  mutate(V4=gsub("chr","",V4),
         liftedId=V1 %&% ":" %&% V2) %>%
  rename(preliftId="V4") %>%
  select(preliftId,liftedId)

head(lift_map$preliftId)

annotation<-fread(input,skip=40,header=T)
str(annotation)
inner_join(lift_map,annotation,by=c('preliftId'='#Uploaded_variation')) %>%
  select(-preliftId) %>%
  rename(`#Uploaded_variation`="liftedId") %>%
  mutate(Location=`#Uploaded_variation`) %>%
  fwrite(.,output,sep='\t')



