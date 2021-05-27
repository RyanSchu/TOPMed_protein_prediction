##################################################################################################

#SETUP ENVIRONMENT
cat("SETUP ENVIRONMENT\n")
##################################################################################################

library(dplyr)
library(qvalue)
library(data.table)
library(ggplot2) 


##################################################################################################

#DEFINE FUNCTIONS
cat("DEFINE FUNCTIONS\n")
##################################################################################################


"%&%" = function(a,b) paste(a,b,sep="")


get_gene_corr<-function(gene_name,m1,m2)
{
  #gene name is the name of the gene
  #colmatrix is a matrix containing genes as columns
  #rowmatrix is a matrix containing genes as rows
  #skipcol = The first n columns of the row matrix do not contain values, skip these for correlation test  

  predicted_exp<-m1 %>% 
    select(gene_name) %>% ##select the gene
    unlist() %>% unname()
  
  measured_exp<-m2 %>% 
    select(gene_name) %>% 
    unlist() %>% 
    unname()
  
  correlation<-cor.test(measured_exp,predicted_exp, method = "spearman")
  
  expression_corr<-list()
  expression_corr[["gene_id"]]<-gene_name
  expression_corr[["estimate"]]<-correlation$estimate
  expression_corr[["p.value"]]<-correlation$p.value
  
  return(expression_corr)
}


##################################################################################################

#SET GLOBAL VARIABLES
cat("SET GLOBAL VARIABLES\n")
##################################################################################################
predicted_dir<-"/home/rschubert1/scratch/TOPMed_Proteome/10INTERVAL_prediction/valid_models/"
file_list<-list.files(path=predicted_dir,pattern="*predicted_expression.txt")

##################################################################################################

#READ IN & PROCESS DATA
cat("READ IN & PROCESS DATA\n")
##################################################################################################

observed_expression<-fread("zcat /data/rschubert1/topmed_proteome_ashg_prep/protein_expr/renamed_protein_expr_raw.txt.gz", header = T, sep = '\t',stringsAsFactors = F) %>%
  rename_at(vars(1),function(x){return("gene_id")})
obsGenes<-gsub("\\.[0-9]+","",colnames(observed_expression))
colnames(observed_expression)<-obsGenes
print(head(obsGenes))
for (i in 1:length(file_list))  {#How well does each model replicate
    predicted_expression<-fread(predicted_dir %&% file_list[i],header=T,stringsAsFactors = F)
    if( sum(predicted_expression$IID == observed_expression$gene_id) != length(sum(predicted_expression$IID)))
    {
      cat("samples not in matching order. Sorting predicted expression to match observed_expression$gene_id.\n")
      predicted_expression<-predicted_expression[match(observed_expression$gene_id,predicted_expression$IID),]
    }
    
    predGenes<-gsub("\\.[0-9]+","",colnames(predicted_expression))
    colnames(predicted_expression)<-predGenes
    print(head(predGenes))
    gene_list<-data.frame(gene_id=base::intersect(obsGenes,predGenes),stringsAsFactors = F)
    
    cat("Existing genes: ",dim(gene_list)[1],"\n")
    if (dim(gene_list)[1] <= 1) { 
      print("No intersection between predicted and observed sets. Make sure ids are formatted properly and matrices are the proper orientation"); 
      next}
    
    ##################################################################################################
    
    #ANALYZE DATA
    cat("ANALYZE DATA\n")
    ##################################################################################################
    # str(predicted_expression)
    # str(observed_expression)
    predictive_correlations<-sapply(X=gene_list$gene_id,FUN=get_gene_corr,m1=predicted_expression,m2=observed_expression,simplify=T,USE.NAMES = T)
    # str(predictive_correlations)
    predictive_correlations<-data.frame(gene_id=unlist(predictive_correlations[1,]),estimate=unlist(predictive_correlations[2,]),p.value=unlist(predictive_correlations[3,]))
    
    # str(predictive_correlations)
    
    ##################################################################################################
    
    #WRITE OUT CORRELATION DATA
    cat("WRITE OUT CORRELATION DATA\n")
    ##################################################################################################
    out_tag<-gsub("predicted_expression.txt","",file_list[i])
    fwrite(predictive_correlations,"/home/rschubert1/scratch/TOPMed_Proteome/10INTERVAL_prediction/valid_model_correlations/" %&% out_tag %&% "correlation.txt",col.names = T,sep='\t')
}#close

