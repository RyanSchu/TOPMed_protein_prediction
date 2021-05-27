#translate id names
library(dplyr)
library(data.table)

mapping<-fread("/data/rschubert1/topmed_proteome_ashg_prep/protein_anno_EGAF00001998175/joint_id_to_somamer_map.tsv")
protein_expr<-fread("zcat /data/rschubert1/topmed_proteome_ashg_prep/protein_expr/protein_expr_raw.txt.gz")

original_names<-data.frame(original=colnames(protein_expr)[-1])
dictionary<-mapping %>% 
  select_at(vars(c(1,10)))
dictionary<-dictionary[!(duplicated(dictionary$joint_id) | duplicated(dictionary$joint_id, fromLast = TRUE)), ]
dictionary<-dictionary[!(duplicated(dictionary$SOMAMER_ID) | duplicated(dictionary$SOMAMER_ID, fromLast = TRUE)), ]
translation<-inner_join(original_names,dictionary,by=c("original"="SOMAMER_ID"))



old_names<-translation %>% select_at(vars(1)) %>% unlist() %>% unname()
old_names<-c(colnames(protein_expr)[1],old_names)
protein_expr<- protein_expr %>% select(one_of(old_names))
#str(protein_expr)
translation<-translation %>% filter(original %in% colnames(protein_expr))

new_names<-translation %>% select_at(vars(2)) %>% unlist() %>% unname()
new_names<-c(colnames(protein_expr)[1],new_names)
str(new_names)

colnames(protein_expr)<-new_names
fwrite(protein_expr,"/data/rschubert1/topmed_proteome_ashg_prep/protein_expr/renamed_protein_expr_raw.txt",sep='\t')