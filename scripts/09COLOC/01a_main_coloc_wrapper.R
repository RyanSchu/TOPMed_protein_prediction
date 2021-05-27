# Authors: Ryan Schubert and Elyse Geoffroy
# This script runs the bse, p, and tstat methods of coloc for QTL and GWAS summary statistics. 
# Performing coloc with a LD matrix for the bse method is an option. 

suppressPackageStartupMessages(library(coloc))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))
"%&%" = function(a,b) paste(a,b,sep="")

check.fread<-function(file,header=T,stringsAsFactors=F){
  #Author Ryan Schubert
  #Meant to be run in unix environments with the zcat program installed
  #Not for the R environ, but useful for the client
  if (grepl("\\.gz$",file)){
    data<-tryCatch({
      fread("zcat " %&% file,header=header,stringsAsFactors = stringsAsFactors)
      },error=function(cond){
        message("returned error") 
        mesasge(cond)
        return(NA)
        })
    if (is.na(data)){
      data<-fread(file,header=header,stringsAsFactors = stringsAsFactors)
      return(data)
    } else {
      return(data)
    }                   
  } else {
    fread(file,header=header,stringsAsFactors = stringsAsFactors)
  }
}


get_gene_list<-function(eqtl_df,geneCol=1){
  #Author Ryan Schubert
  #takes in the data frame of eqtl summary stats and
  #extracts the column with the gene names in them,
  #returns as a vector
  eqtl_df[,geneCol] %>% unlist %>% unname %>% unique()
}


get_eqtl_snps<-function(gene_id,eqtl_df,snpCol=2,geneCol=1){
  #Author Ryan Schubert
  #take in the eqtl data frame
  #returns the snps in the data frame
  #If the gene is not found in the eqtl df then return an empty list
  #This is a separate function from the effects and pval function because the maf is not always included in the eqtl data, though it is common practice
  eqtl_df<- eqtl_df[eqtl_df[,geneCol] == gene_id,]
  # print(gene_id)
  # print(eqtl_df)
  if (nrow(eqtl_df) == 0) { return(list())}
  snp<-eqtl_df[,snpCol] %>% unlist %>% unname
  # print(snp)
  return(list(snp=snp))
}

get_gene_eqtl_effects<-function(gene_id,eqtl_df,betaCol,seCol,geneCol=1,snpCol,snpList){
  #Author Ryan Schubert, Edits by Elyse Geoffroy
  #take in the eqtl data frame
  #filters to the gene in question
  #returns the effect sizes and variances for each snp
  #If the gene is not found in the eqtl df then return an empty list
  print(gene_id)
  eqtl_df<-eqtl_df[eqtl_df[,geneCol] == gene_id & eqtl_df[,snpCol] %in% snpList,]
  if (nrow(eqtl_df) == 0) { return(list())}
  beta <- eqtl_df[,betaCol] %>% unlist %>% unname
  var <- eqtl_df[,seCol]^2 %>% unlist %>% unname
  snp <- eqtl_df[,snpCol] %>% unlist %>% unname
  return(list(beta=beta,var=var,snp=snp))
}

get_gene_eqtl_pvalue<-function(gene_id,eqtl_df,pvalCol,geneCol=1,snpCol,snpList){
  #Author Ryan Schubert, Edits by Elyse Geoffroy
  #take in the eqtl data frame
  #filters to the gene in question
  #returns the pvalue of each snp
  #If the gene is not found in the eqtl df then return an empty list
  #separate function from the effects function because for coloc purposes you only need one or the other, not both
  eqtl_df<-eqtl_df[eqtl_df[,snpCol] %in% snpList,]
  eqtl_df <- eqtl_df[eqtl_df[,geneCol]==gene_id, ]
  if (nrow(eqtl_df) == 0) { return(list())}
  p<-eqtl_df[,..pvalCol] %>% unlist %>% unname
  return(list(p=p))
}

get_gene_eqtl_tstat<-function(gene_id,eqtl_df,tstatCol,geneCol=1,snpCol,snpList){
  #Author Ryan Schubert
  #take in the eqtl data frame
  #filters to the gene in question
  #returns the pvalue of each snp
  #If the gene is not found in the eqtl df then return an empty list
  #separate function from the effects function because for coloc purposes you only need one or the other, not both
  eqtl_df<-eqtl_df[eqtl_df[,geneCol] == gene_id & eqtl_df[,snpCol] %in% snpList,]
  if (nrow(eqtl_df) == 0) { return(list())}
  t<-eqtl_df[,tstatCol] %>% unlist %>% unname
  return(list(t=t))
}

get_gene_eqtl_maf<-function(gene_id,eqtl_df,mafCol=3,geneCol=1,snpCol,snpList, sigsnps=NULL){
  #Author Ryan Schubert, Edits by Elyse Geoffroy
  #take in the eqtl data frame
  #filters to the gene in question
  #returns the maf of each snp
  #If the gene is not found in the eqtl df then return an empty list
  #This is a separate function from the effects and pval function because the maf is not always included in the eqtl data, though it is common practice
  eqtl_df<-eqtl_df[eqtl_df[,geneCol] == gene_id & eqtl_df[,snpCol] %in% snpList,]
  if(!is.null(sigsnps)){
        eqtl_df <- eqtl_df[eqtl_df[,snpCol] %in% sigsnps, ]
	# print(head(eqtl_df))
  }
  if (nrow(eqtl_df) == 0) { return(list())}
  maf<-eqtl_df[,mafCol] %>% unlist %>% unname
  return(list(maf=maf))
}


get_gwas_snps<-function(gwas_df,snpCol=1){
  #Author Ryan Schubert
  #take in the gwas data frame
  #returns the snps in the data frame
  #If the gene is not found in the gwas df then return an empty list
  #This is a separate function from the effects and pval function because the maf is not always included in the gwas data, though it is common practice
  snp<-gwas_df[,snpCol] %>% unlist %>% unname
  return(list(snp=snp))
}

get_gwas_effects<-function(gwas_df,betaCol,seCol,snpCol,snpList){
  #Author Ryan Schubert
  #take in the gwas data frame
  #returns the effect sizes and variances for each snp
  #If the gene is not found in the gwas df then return an empty list
  #print(snpList)
  gwas_df<-gwas_df[gwas_df[,snpCol] %in% snpList,]
  if (nrow(gwas_df) == 0) { return(list())}
  beta<-gwas_df[,betaCol] %>% unlist %>% unname
  var<-gwas_df[,seCol]^2 %>% unlist %>% unname
  snp<-gwas_df[,snpCol] %>% unlist %>% unname
  return(list(beta=beta,var=var,snp=snp))
}

get_gwas_pvalue<-function(gwas_df,pvalCol,snpCol,snpList){
  #Author Ryan Schubert
  #take in the gwas data frame
  #returns the pvalue of each snp
  #If the gene is not found in the gwas df then return an empty list
  #separate function from the effects function because for coloc purposes you only need one or the other, not both
  gwas_df<-gwas_df[gwas_df[,snpCol] %in% snpList,]
  if (nrow(gwas_df) == 0) { return(list())}
  p<-gwas_df[,pvalCol] %>% unlist %>% unname
  return(list(p=p))
}

get_gwas_tstat<-function(gwas_df,tstatCol,snpCol,snpList){
  #Author Ryan Schubert
  #take in the gwas data frame
  #returns the pvalue of each snp
  #If the gene is not found in the gwas df then return an empty list
  #separate function from the effects function because for coloc purposes you only need one or the other, not both
  gwas_df<-gwas_df[gwas_df[,snpCol] %in% snpList,]
  if (nrow(gwas_df) == 0) { return(list())}
  t<-gwas_df[,tstatCol] %>% unlist %>% unname
  return(list(t=t))
}

get_gwas_maf<-function(gwas_df,mafCol,snpCol,snpList, sigsnps=NULL){
  #Author Ryan Schubert, Edits by Elyse Geoffroy
  #take in the gwas data frame
  #returns the maf of each snp
  #If the gene is not found in the gwas df then return an empty list
  #This is a separate function from the effects and pval function because the maf is not always included in the gwas data, though it is common practice
  gwas_df<-gwas_df[gwas_df[,snpCol] %in% snpList,]
  if(!is.null(sigsnps)){
	gwas_df <- gwas_df[gwas_df[,snpCol] %in% sigsnps, ]
  }
  if (nrow(gwas_df) == 0) { return(list())}
  maf<-gwas_df[,..mafCol] %>% unlist %>% unname
  return(list(maf=maf))
}

return_ordered_snps<-function(snpList,eqtl_df,snpCol){
  #Author Ryan Schubert
  #Take in the list of snps of interest
  #return a list of snps that are present in the eqtl df in the order that they appear in the eqtl df
  snps<-eqtl_df[eqtl_df[,snpCol] %in% snpList,snpCol] %>% unlist %>% unname
  return(list(snps=snps))
}


main<-function(eqtl,gwas,mode="bse", gene_list=NULL, directory='',
               eqtlGeneCol=NULL,eqtlSNPCol=NULL,eqtlMAFCol=NULL,eqtlPvalCol=NULL,eqtlBetaCol=NULL,eqtlSeCol=NULL, eqtlTstatCol=NULL,eqtlSampleSize=NULL,
               gwasSNPCol=NULL,gwasMAFCol=NULL,gwasPvalCol=NULL,gwasBetaCol=NULL,gwasSeCol=NULL, gwasTstatCol=NULL,gwasSampleSize=NULL,rule="H4>0.5",
               outFile=NULL, LD=NULL, ld_header=NULL,
               method=NULL,p12=1e-5,p1=1e-4,p2=1e-4,pthr=1e-6,r2thr=0.01,maxhits=3,cmode="iterative"){
  if(is.null(eqtlSampleSize) || is.null(gwasSampleSize)) {print("Sample sizes not set. Please make sure both GWAS and QTL sample sizes are set. Exiting."); return()}
  if( mode  == "bse") {
    cat("Running coloc using the beta/standard error method\n")
    if(is.null(eqtlBetaCol) || is.null(eqtlSeCol)) {print("eQTL Beta/Se column not set. Please set these and try again. Exiting"); return() }
    if(is.null(gwasBetaCol) || is.null(gwasSeCol)) {print("GWAS Beta/Se column not set. Please set these and try again. Exiting"); return() }
  } else if ( mode == "p" ){
    cat("Running coloc using the pvalue method\n")
    if(is.null(eqtlPvalCol)) {print("eQTL Pval column not set. Please set this and try again. Exiting"); return() }
    if(is.null(gwasPvalCol)) {print("GWAS Pval column not set. Please set this and try again. Exiting"); return() }
  } else if ( mode == "t" ){
      cat("Running coloc using the t-stat method\n")
      if(is.null(eqtlTstatCol)) {print("eQTL Tstat column not set. Please set this and try again. Exiting"); return() }
      if(is.null(gwasTstatCol)) {print("GWAS Tstat column not set. Please set this and try again. Exiting"); return() }
  } else {
    cat("Mode not recognized. Please choose one of [ bse, p, t ]. Exiting\n"); return()
  }

  #uncomment the follwing lines when running in a Unix/linux environ
   cat("Reading in data\n")
   eqtldf<-as.data.frame(check.fread(eqtl,header=T,stringsAsFactors = F)) %>% distinct()
   gwasdf<-as.data.frame(check.fread(gwas,header=T,stringsAsFactors = F)) %>% distinct()

  #comment out the following lines when running in Unix/Linux environ
#  cat("Reading in data\n")
  if(nrow(eqtldf) == 0){
  	eqtldf<-read.table(eqtl,header=T,stringsAsFactors = F)
  	eqtldf <- as.data.frame(unique(eqtldf))
  }
  #print(head(eqtldf))
  if(nrow(gwasdf) == 0){
	gwasdf<-read.table(gwas,header=T,stringsAsFactors = F)
        gwasdf <- as.data.frame(unique(gwasdf))
  }
  #print(head(gwasdf))

#  cat("Getting gene list\n")
  #gene_list<-get_gene_list(eqtldf)
  ngenes<-length(gene_list)
  if ( ngenes == 0 ) {print("0 genes found. Exiting.\n"); return()}
#  cat(ngenes,"genes found\n")
  if ( !is.null(outFile) ){
    names<-list(hit2=character(),hit1=character(),nsnps=numeric(),PP.H0.abf=numeric(),PP.H1.abf=numeric(),PP.H2.abf=numeric(),PP.H3.abf=numeric(),PP.H4.abf=numeric(),best1=character(),best2=character(),best4=character(),hit1.margz=numeric(),hit2.margz=numeric(),gene=character())
    write.table(names,outFile,sep='\t',quote=F)
    #names2 <- list(gene_id=character(), p0=numeric(),p1=numeric(),p2=numeric(),p3=numeric(),p4=numeric())
    #write.table(names2, paste('gene_formatted_', outFile, sep = ''), sep='\t',quote=F)
  }

  if ( mode == "bse") {
    cat("Processing GWAS\n")

    GWAS_snps<-get_gwas_snps(gwas_df = gwasdf,snpCol=gwasSNPCol)
    for (i in 1:ngenes){
      gene<-gene_list[i]
      cat("Processing gene",gene," ",i,"/",ngenes,"\n")
      QTL_snps<-get_eqtl_snps(gene_id=gene,eqtl_df=eqtldf,snpCol=eqtlSNPCol,geneCol=eqtlGeneCol)
      #print(QTL_snps)
      intersection<-base::intersect(GWAS_snps$snp,QTL_snps$snp)
      nsnps<-length(intersection)
      if (nsnps == 0) { print("0 snps present in intersection for this gene. skipping"); next}


      cat(nsnps, "snps found in GWAS and eQTL intersection for this gene\n")
      print(head(intersection))
      GWAS_effects<-get_gwas_effects(gwas_df = gwasdf,betaCol = gwasBetaCol,seCol=gwasSeCol,snpCol=gwasSNPCol,snpList=intersection)
      QTL_effects<-get_gene_eqtl_effects(gene_id=unique(gene),eqtl_df = eqtldf,betaCol=eqtlBetaCol,seCol=eqtlSeCol,snpCol=eqtlSNPCol,snpList=intersection)
      print(length(GWAS_effects$beta))
      #print(head(QTL_effects))
      print(length(QTL_effects$beta))
      # fwrite(QTL_effects,"/home/rschubert1/scratch/test_coloc_for_elyse/QTL_effects.txt",sep='\t')
      # fwrite(GWAS_effects,"/home/rschubert1/scratch/test_coloc_for_elyse/GWAS_effects.txt",sep='\t')
      if (length(GWAS_effects$beta) != length(QTL_effects$beta)) { print(gene %&% " Lists of effect sizes of differing lengths. Duplicates SNPs may be present. Please resolve and rerun. Exiting."); next()}
      print('after effect sizes check')

	 # will have to filter out snps from GWAS and eQTL that are in range of ld matrix
        
	if(!is.null(LD)){
        #ld_matrix <- fread(directory %&% gene %&% ', header = F, stringsAsFactors=F)
	if(!is.null(list.files(directory, pattern = paste(gene, '_1Mb_LD.ld.gz', sep =''), full.names=T)[[1]][1])){
        	ld_matrix1 <- list.files(directory, pattern = paste(gene, '_1Mb_LD.ld.gz', sep =''), full.names=T)[[1]][1]
        	ld_matrix <- fread(ld_matrix1, header = F, stringsAsFactors=F)
        	ld_matrix <- as.matrix(ld_matrix)
        	ld_header <- str_replace(ld_matrix1, '.ld.gz', '.snplist')
		ld_header <- fread(ld_header, header = F, sep = ':', stringsAsFactors=F)
		      if (grepl("chr",intersection[1])){
		        ld_header<-mutate(ld_header,V1=if_else(!grepl("chr",V1),"chr" %&% V1,as.character(V1)))
		      }  else {
		        ld_header<-mutate(ld_header,V1=if_else(grepl("chr",V1),gsub("chr","",V1),as.character(V1)))
		      }
        	ld_header <- paste(ld_header$V1, ld_header$V2, sep =':')
        	colnames(ld_matrix) <- as.list(ld_header)
        	rownames(ld_matrix) <- as.list(ld_header)

        	#print(head(ld_matrix))
        	print(dim(ld_matrix))
	} else {cat(gene, 'gene skipped due to lack of ld matrix file when it was requested to use new version of coloc'); next}
	}

      if (!is.null(gwasMAFCol)){
        print("using maf from gwas data set")
        
        maf<-get_gwas_maf(gwas_df=gwasdf,mafCol=gwasMAFCol,snpCol = gwasSNPCol,snpList=intersection, sigsnps = ld_header)
      } else if(!is.null(eqtlMAFCol)){
        print("using maf from QTL data set")
        maf<-get_gene_eqtl_maf(gene_id=gene,eqtl_df = eqtldf,mafCol=eqtlMAFCol,geneCol=eqtlGeneCol,snpCol=eqtlSNPCol,snpList=intersection, sigsnps = ld_header)
      }
      

	
      print("here")
      GWAS_effects <- data.frame(GWAS_effects,stringsAsFactors = F)
      QTL_effects <- data.frame(QTL_effects,stringsAsFactors = F)
      print(str(ld_header))
      print(GWAS_effects$snp[1])
      # if(grepl("chr",ld_header[1])){
      #   GWAS_effects<-mutate(GWAS_effects,snp=if_else(!grepl("chr",snp),"chr" %&% snp,snp))
      #   QTL_effects<-mutate(QTL_effects,snp=if_else(!grepl("chr",snp),"chr" %&% snp,snp))
      # } else {
      #   GWAS_effects<-mutate(GWAS_effects,snp=if_else(grepl("chr",snp),gsub("chr","",snp),snp))
      #   QTL_effects<-mutate(QTL_effects,snp=if_else(grepl("chr",snp),gsub("chr","",snp),snp))
      # }	
	GWAS_effects <- data.frame(GWAS_effects)
	GWAS_effects <- GWAS_effects[GWAS_effects$snp %in% ld_header, ] #GWAS_effects[GWAS_effects[,gwasSNPCol] %in% ld_header, ]
	GWAS_effects$snp <- as.character(GWAS_effects$snp)
	print(nrow(GWAS_effects))

	QTL_effects <- data.frame(QTL_effects)
	QTL_effects <- QTL_effects[QTL_effects$snp %in% ld_header, ] #QTL_effects[QTL_effects[,eqtlSNPCol] %in% sig_snps, ] 
	QTL_effects$snp <- as.character(QTL_effects$snp)
        print(nrow(QTL_effects))

      if(nrow(GWAS_effects) > 0){
      cat('running coloc\n')
      str(list(beta=GWAS_effects$beta, varbeta=GWAS_effects$var,snp=GWAS_effects$snp, N=gwasSampleSize,type="quant"))
      str(list(beta=QTL_effects$beta, varbeta=QTL_effects$var, snp=QTL_effects$snp, N=eqtlSampleSize,type="quant"))
      
      coloc_result<-coloc.signals(dataset1=list(beta=GWAS_effects$beta, varbeta=GWAS_effects$var,snp=GWAS_effects$snp, N=gwasSampleSize,type="quant"),
                    dataset2=list(beta=QTL_effects$beta, varbeta=QTL_effects$var, snp=QTL_effects$snp, N=eqtlSampleSize,type="quant"),
                    MAF=maf$maf, LD=ld_matrix,
                    method=method,p12=p12,p1=p1,p2=p2,pthr=pthr,r2thr=r2thr,mode=cmode,maxhits=maxhits)
      print("coloc result finished.")
      # if ( coloc_result$summary[6]> 0.5 ){
      if (nrow(coloc_result$summary)>0){
          pdfname<-as.character(unlist(strsplit(outFile, split = ".txt"))) #was split by '_'
          pdf(paste0(pdfname[1],"-t-",gene,".pdf")) #was 2
	  sensitivity(coloc_result,rule=rule)
          dev.off()
      }
      summary<-coloc_result$summary
      summary$gene<-gene
      }
      
      if ( !is.null(outFile) ){
        write.table(summary,outFile,append=T,sep='\t',col.names=F,quote=F,row.names=F)
      }
#      print(summary)
    }
 } else if (mode == "p"){
    cat("Processing GWAS\n")

    GWAS_snps<-get_gwas_snps(gwas_df = gwasdf,snpCol=gwasSNPCol)

    for (i in 1:ngenes){
      gene<-gene_list[i]
      cat("Processing gene",gene," ",i,"/",ngenes,"\n")
      QTL_snps<-get_eqtl_snps(gene_id=gene,eqtl_df=eqtldf,snpCol=eqtlSNPCol,geneCol=eqtlGeneCol)
      intersection<-base::intersect(GWAS_snps$snp,QTL_snps$snp)
      nsnps<-length(intersection)
      if (nsnps == 0) { print("0 snps present in intersection for this gene. skipping"); next}
      cat(nsnps, "snps found in GWAS and eQTL intersection for this gene\n")
      GWAS_effects<-get_gwas_pvalue(gwas_df = gwasdf, pvalCol = gwasPvalCol,snpCol=gwasSNPCol,snpList=intersection)
      QTL_effects<-get_gene_eqtl_pvalue(gene_id=gene,eqtl_df = eqtldf, pvalCol=eqtlPvalCol,snpCol=eqtlSNPCol,snpList=intersection)
      if (length(GWAS_effects$p) != length(QTL_effects$p)) { print("List of effect size of differing lengths. Duplicates SNPs may be present. Please resolve and rerun. Exiting."); return()}
	

      if(!is.null(LD)){
        #ld_matrix <- fread(directory %&% gene %&% ', header = F, stringsAsFactors=F)
      	if(!is.null(list.files(directory, pattern = paste(gene, '_1Mb_LD.ld.gz', sep =''), full.names=T)[[1]][1])){
        	ld_matrix1 <- list.files(directory, pattern = paste(gene, '_1Mb_LD.ld.gz', sep =''), full.names=T)[[1]][1]
        	print(ld_matrix1)
        	ld_matrix <- fread(ld_matrix1, header = F, stringsAsFactors=F)
        	ld_matrix <- as.matrix(ld_matrix)
        	ld_header <- str_replace(ld_matrix1, '.ld.gz', '.snplist')
		      ld_header <- fread(ld_header, header = F, sep = ':', stringsAsFactors=F)
        	ld_header <- paste(ld_header$V1, ld_header$V2, sep =':')
        	colnames(ld_matrix) <- as.list(ld_header)
        	rownames(ld_matrix) <- as.list(ld_header)

        	#print(head(ld_matrix))
        	print(dim(ld_matrix))
	        } else {cat(gene, 'gene skipped due to lack of ld matrix file when it was requested to use new version of coloc')}
         }

      if (!is.null(gwasMAFCol)){
        print("using maf from gwas data set")
        maf<-get_gwas_maf(gwas_df=gwasdf,mafCol=gwasMAFCol,snpCol = gwasSNPCol,snpList=intersection)
      } else if(!is.null(eqtlMAFCol)){
        print("using maf from QTL data set")
        print("here")
        maf<-get_gene_eqtl_maf(gene_id=gene,eqtl_df = eqtldf,mafcol=eqtlMAFCol,geneCol=eqtlGeneCol,snpCol=eqtlSNPCol,snpList=intersection)
      }
      # str(maf$maf); str(GWAS_effects$p); str(QTL_effects$p)
         print("here")
	      GWAS_effects <- data.frame(GWAS_effects)
	      QTL_effects <- data.frame(QTL_effects)
	      print(str(ld_header))
	      print(GWAS_effects$snp[1])
	      if(grepl("chr",ld_header[1])){
	        GWAS_effects<-mutate(GWAS_effects,snp=if_else(!grepl("chr",snp),"chr" %&% snp,snp))
	        QTL_effects<-mutate(QTL_effects,snp=if_else(!grepl("chr",snp),"chr" %&% snp,snp))
	      } else {
	        GWAS_effects<-mutate(GWAS_effects,snp=if_else(grepl("chr",snp),gsub("chr","",snp),snp))
	        QTL_effects<-mutate(QTL_effects,snp=if_else(grepl("chr",snp),gsub("chr","",snp),snp))
	      }
        GWAS_effects <- GWAS_effects[GWAS_effects$snp %in% ld_header, ] #GWAS_effects[GWAS_effects[,gwasSNPCol] %in% ld_header, ]
        GWAS_effects$snp <- as.character(GWAS_effects$snp)
        print(nrow(GWAS_effects))
         print("now here")
        
        QTL_effects <- QTL_effects[QTL_effects$snp %in% ld_header, ] #QTL_effects[QTL_effects[,eqtlSNPCol] %in% sig_snps, ]
        QTL_effects$snp <- as.character(QTL_effects$snp)
        print(nrow(QTL_effects))
         print("here next")
      coloc_result<-coloc.signals(dataset1=list(pvalues=GWAS_effects$p, N=gwasSampleSize,type="quant"),
                                  dataset2=list(pvalues=GWAS_effects$p, N=gwasSampleSize,type="quant"),
                                  MAF=maf$maf, LD=ld_matrix,
                                  method=method,p12=p12,p1=p1,p2=p2,pthr=pthr,r2thr=r2thr,mode=cmode,maxhits=maxhits)
        print("finally here")
      # if ( coloc_result$summary[6]> 0.5 ){
      if (nrow(coloc_result$summary)>0){
        pdfname<-as.character(unlist(strsplit(outFile, split = "_")))
	pdf(paste0(pdfname[2],"-t-",gene,".pdf"))
        sensitivity(coloc_result,rule=rule)
        dev.off()
      }
      summary<-coloc_result$summary
      summary$gene<-gene
      # summary<-bind_cols(summary)
      # str(summary)
      if ( !is.null(outFile) ){
        write.table(summary,outFile,append=T,sep='\t',col.names=F,quote=F,row.names=F)
      }
    }
      #      print(summary)
    } else if (mode == "t"){
      cat("Processing GWAS\n")

      GWAS_snps<-get_gwas_snps(gwas_df = gwasdf,snpCol=gwasSNPCol)

      for (i in 1:ngenes){
        gene<-gene_list[i]
        cat("Processing gene",gene," ",i,"/",ngenes,"\n")
        QTL_snps<-get_eqtl_snps(gene_id=gene,eqtl_df=eqtldf,snpCol=eqtlSNPCol,geneCol=eqtlGeneCol)
        intersection<-base::intersect(GWAS_snps$snp,QTL_snps$snp)
        nsnps<-length(intersection)
        if (nsnps == 0) { print("0 snps present in intersection for this gene. skipping"); next}
        cat(nsnps, "snps found in GWAS and eQTL intersection for this gene\n")
        GWAS_effects<-get_gwas_tstat(gwas_df = gwasdf, tstatCol = gwasTstatCol,snpCol=gwasSNPCol,snpList=intersection)
        QTL_effects<-get_gene_eqtl_tstat(gene_id=gene,eqtl_df = eqtldf, tstatCol=eqtlTstatCol,snpCol=eqtlSNPCol,snpList=intersection)
        if (length(GWAS_effects$t) != length(QTL_effects$t)) { print("List of effect size of differing lengths. Duplicates SNPs may be present. Please resolve and rerun. Exiting."); return()}
        if(!is.null(LD)){
        #ld_matrix <- fread(directory %&% gene %&% ', header = F, stringsAsFactors=F)
	if(!is.null(list.files(directory, pattern = paste(gene, '_1Mb_LD.ld.gz', sep =''), full.names=T)[[1]][1])){
        	ld_matrix1 <- list.files(directory, pattern = paste(gene, '_1Mb_LD.ld.gz', sep =''), full.names=T)[[1]][1]
        	ld_matrix <- fread(ld_matrix1, header = F, stringsAsFactors=F)
        	ld_matrix <- as.matrix(ld_matrix)
        	ld_header <- str_replace(ld_matrix1, '.ld.gz', '.snplist')
		ld_header <- fread(ld_header, header = F, sep = ':', stringsAsFactors=F)
        	ld_header <- paste(ld_header$V1, ld_header$V2, sep =':')
        	colnames(ld_matrix) <- as.list(ld_header)
        	rownames(ld_matrix) <- as.list(ld_header)

        	#print(head(ld_matrix))
        	print(dim(ld_matrix))
	} else {cat(gene, 'gene skipped due to lack of ld matrix file when it was requested to use new version of coloc')}
	}

        if (!is.null(gwasMAFCol)){
          print("using maf from gwas data set")
          maf<-get_gwas_maf(gwas_df=gwasdf,mafCol=gwasMAFCol,snpCol = gwasSNPCol,snpList=intersection)
        } else if(!is.null(eqtlMAFCol)){
          print("using maf from QTL data set")
          maf<-get_gene_eqtl_maf(gene_id=gene,eqtl_df = eqtldf,mafcol=eqtlMAFCol,geneCol=eqtlGeneCol,snpCol=eqtlSNPCol,snpList=intersection)
        }
        GWAS_effects <- data.frame(GWAS_effects)
        GWAS_effects <- GWAS_effects[GWAS_effects$snp %in% ld_header, ] #GWAS_effects[GWAS_effects[,gwasSNPCol] %in% ld_header, ]
        GWAS_effects$snp <- as.character(GWAS_effects$snp)
        print(nrow(GWAS_effects))

        QTL_effects <- data.frame(QTL_effects)
        QTL_effects <- QTL_effects[QTL_effects$snp %in% ld_header, ] #QTL_effects[QTL_effects[,eqtlSNPCol] %in% sig_snps, ]
        QTL_effects$snp <- as.character(QTL_effects$snp)
        print(nrow(QTL_effects))

        coloc_result <- coloc.abf(dataset1=list(pval=GWAS_effects$t, N=gwasSampleSize,type="quant"),
                                  dataset2=list(pval=QTL_effects$t, N=eqtlSampleSize,type="quant"),
                                  MAF=maf$maf, LD=ld_matrix)
        if (coloc_result$summary[6]>0.5){
        if (coloc_result$summary[1]>1){
          pdfname<-as.character(unlist(strsplit(outFile, split = "_")))
          pdf(paste0(pdfname[2],"-t-",gene,".pdf"))
          sensitivity(coloc_result,rule=rule)
          dev.off()
        }
      }
        summary<-coloc_result$summary
        summary$gene<-gene
        summary<-bind_cols(summary)
        # str(summary)
        if ( !is.null(outFile) ){
          write.table(summary,outFile,append=T,sep='\t',col.names=F,quote=F,row.names=F)
        }
        #      print(summary)
      }
    }
}
