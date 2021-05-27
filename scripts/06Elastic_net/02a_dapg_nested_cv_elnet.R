suppressMessages(library(dplyr))
suppressMessages(library(glmnet))
suppressMessages((library(reshape2)))
suppressMessages(library(methods))
suppressMessages(library(doMC))
suppressMessages(library(doRNG))
suppressMessages(library(tidyr))
"%&%" <- function(a,b) paste(a,b, sep = "")
#signal #

get_filtered_snp_annot <- function(snp_annot_file_name) {
  snp_annot <- read.table(snp_annot_file_name, header = T, stringsAsFactors = F) %>%
    filter(!((refAllele == 'A' & effectAllele == 'T') |
               (refAllele == 'T' & effectAllele == 'A') |
               (refAllele == 'C' & effectAllele == 'G') |
               (refAllele == 'G' & effectAllele == 'C')) &
             !(is.na(rsid))) %>%
    distinct(varID, .keep_all = TRUE)
  snp_annot
}


get_maf_filtered_genotype <- function(genotype_file_name,  maf, samples,snp_annot) {
  gt_df <- read.table(genotype_file_name, header = T, stringsAsFactors = F, row.names = 1)
  gt_df <- gt_df[rownames(gt_df) %in% snp_annot$rsid,]
  gt_df <- gt_df[,samples] %>% t() %>% as.data.frame()
  effect_allele_freqs <- colMeans(gt_df) / 2
  gt_df <- gt_df[,which((effect_allele_freqs >= maf) & (effect_allele_freqs <= 1 - maf))]
  gt_df
}

get_gene_annotation <- function(gene_annot_file_name, chrom, gene_types=c('protein_coding', 'pseudogene', 'lincRNA','aptamer')){
  gene_df <- read.table(gene_annot_file_name, header = TRUE, stringsAsFactors = FALSE) %>%
    filter((chr == chrom) & gene_type %in% gene_types)
  gene_df
}

get_gene_type <- function(gene_annot, gene) {
  filter(gene_annot, gene_id == gene)$gene_type
}

get_gene_expression <- function(gene_expression_file_name, gene_annot) {
  expr_df <- as.data.frame(t(read.table(gene_expression_file_name, header = T, stringsAsFactors = F, row.names = 1)))
  expr_df <- expr_df %>% select(one_of(intersect(gene_annot$gene_id, colnames(expr_df))))
  expr_df
}

get_gene_coords <- function(gene_annot, gene) {
  row <- gene_annot[which(gene_annot$gene_id == gene),]
  c(row$start, row$end)
}

get_cis_genotype <- function(gt_df, snp_annot, coords, cis_window) {
  snp_info <- snp_annot %>% filter((pos >= (coords[1] - cis_window) & !is.na(rsid)) & (pos <= (coords[2] + cis_window)))
  if (nrow(snp_info) == 0)
    return(NA)
  cis_gt <- gt_df %>% select(one_of(intersect(snp_info$varID, colnames(gt_df))))
  column_labels <- colnames(cis_gt)
  row_labels <- rownames(cis_gt)
  # Convert cis_gt to a matrix for glmnet
  cis_gt <- matrix(as.matrix(cis_gt), ncol=ncol(cis_gt)) # R is such a bad language.
  colnames(cis_gt) <- column_labels
  rownames(cis_gt) <- row_labels
  cis_gt
}###This function is now unused in this script

get_covariates <- function(covariate_file_name, samples) {
  cov_df <- tryCatch({read.table(covariate_file_name, header = TRUE, stringsAsFactors = FALSE, row.names = 1)} ,
                     error=function(cond){
                       cond
                     })
  if(inherits(cov_df, "error")) {
    return(NA)
  } else{
    cov_df <- cov_df[,samples] %>% t() %>% as.data.frame()
    cov_df
  }
}###added a try catch error handler for the covariates since our proteome has already been adjusted by covariate. This makes the covariates file optional

generate_fold_ids <- function(n_samples, n_folds=10) {
  n <- ceiling(n_samples / n_folds)
  # cat("fold_ids",n,"\n")
  # cat("nfolds",n_folds,"\n")
  # cat("n_samples",n_samples,"\n")
  fold_ids <- rep(1:n_folds, n)
  sample(fold_ids[1:n_samples])
}

adjust_for_covariates <- function(expression_vec, cov_df) {
  combined_df <- cbind(expression_vec, cov_df)
  expr_resid <- summary(lm(expression_vec ~ ., data=combined_df))$residuals
  expr_resid <- scale(expr_resid, center = TRUE, scale = TRUE)
  expr_resid
}

calc_R2 <- function(y, y_pred) {
  tss <- sum(y**2)
  rss <- sum((y - y_pred)**2)
  1 - rss/tss
}

calc_corr <- function(y, y_pred) {
  sum(y*y_pred) / (sqrt(sum(y**2)) * sqrt(sum(y_pred**2)))
}

read_dapg<-function(dapg_file){ 
  pip_df<-read.table(dapg_file, header = F,stringsAsFactors = F) #%>% separate(V7,into = c("V7","V8"),sep=" ")
  #str(pip_df)
  colnames(pip_df)<-c("rank","snp","PIP","var1","cluster","var2","var3","txt")
  pip_df<-pip_df %>% mutate(ID=gsub(".fm.out","",txt)) %>% select(ID,snp,PIP,cluster)
  pip_df
} ###function just returns all genes in the dapg input, their snps, and their pip

get_dapg_genotype<-function(gt_df,pip_df,gene_name,pip_threshold=0.01,filt_cluster=T) {
   # str(pip_df)
  gene_snp_1 <- pip_df %>% 
    filter(ID==gene_name) %>% 
    filter(PIP > pip_threshold)
  if(filt_cluster==T){
    gene_snp <- gene_snp_1 %>% group_by(cluster) %>% top_n(n=1,wt=PIP) %>% filter(cluster != -1)
  }else {
    gene_snp<-gene_snp_1
  }
   # str(colnames(gt_df))
   # str(gene_snp)
  if (nrow(gene_snp) == 0)
    return(NA)
  dapgenotypes <- gt_df %>% select(one_of(intersect(gene_snp$snp, colnames(gt_df)))) #look I made a pun
  if (ncol(dapgenotypes) == 0)
    return(NA)
  column_labels <- colnames(dapgenotypes)
  row_labels <- rownames(dapgenotypes)
  dapgenotypes <- matrix(as.matrix(dapgenotypes), ncol=ncol(dapgenotypes)) # R is such a bad language.
  colnames(dapgenotypes) <- column_labels
  rownames(dapgenotypes) <- row_labels
  dapgenotypes
}### Processes the dapg dat frame and selects the top snp in each cluster for a gene that meets PIP threshold

get_dapg_penalty<-function(dapgenotypes,pip_df,gene_name){
  obssnps<-colnames(dapgenotypes)
  pip_df %>% filter(ID==gene_name) %>% filter(snp %in% obssnps)
  pip_df<-pip_df[match(obssnps, pip_df$snp),]
  penalty<-pip_df %>% mutate(penalty=1-PIP) %>% select(penalty) %>% unlist()
  penalty
}### Gets the penalty scores from the dapg PIPs. penalty is 1-PIP

nested_cv_elastic_net_perf <- function(x, y, n_samples, n_train_test_folds, n_k_folds, alpha,penalty) { ### added penalty factor 
  # Gets performance estimates for k-fold cross-validated elastic-net models.
  # Splits data into n_train_test_folds disjoint folds, roughly equal in size,
  # and for each fold, calculates a n_k_folds cross-validated elastic net model. Lambda parameter is
  # cross validated. Then get performance measures for how the model predicts on the hold-out
  # fold. Get the coefficient of determination, R^2, and a p-value, where the null hypothesis
  # is there is no correlation between prediction and observed.
  #
  # The mean and standard deviation of R^2 over all folds is then reported, and the p-values
  # are combined using Fisher's method.
  R2_folds <- rep(0, n_train_test_folds)
  corr_folds <- rep(0, n_train_test_folds)
  zscore_folds <- rep(0, n_train_test_folds)
  pval_folds <- rep(0, n_train_test_folds)
  # Outer-loop split into training and test set.
  train_test_fold_ids <- generate_fold_ids(n_samples, n_folds=n_train_test_folds)
  for (test_fold in 1:n_train_test_folds) {
    train_idxs <- which(train_test_fold_ids != test_fold)
    test_idxs <- which(train_test_fold_ids == test_fold)
    x_train <- x[train_idxs, ]
    y_train <- y[train_idxs]
    x_test <- x[test_idxs, ]
    y_test <- y[test_idxs]
    # Inner-loop - split up training set for cross-validation to choose lambda.
    cv_fold_ids <- generate_fold_ids(length(y_train), n_k_folds)
    y_pred <- tryCatch({
      # Fit model with training data.
      fit <- cv.glmnet(x_train, y_train, nfolds = n_k_folds, alpha = alpha, type.measure='mse', foldid = cv_fold_ids, parallel = TRUE,penalty.factor=penalty) ### added penalty.factor=penalty
      # Predict test data using model that had minimal mean-squared error in cross validation.
      predict(fit, x_test, s = 'lambda.min')},
      # if the elastic-net model did not converge, predict the mean of the y_train (same as all non-intercept coef=0)
      error = function(cond) rep(mean(y_train), length(y_test)))
    R2_folds[test_fold] <- calc_R2(y_test, y_pred)
    # Get p-value for correlation test between predicted y and actual y.
    # If there was no model, y_pred will have var=0, so cor.test will yield NA.
    # In that case, give a random number from uniform distribution, which is what would
    # usually happen under the null.
    corr_folds[test_fold] <- ifelse(sd(y_pred) != 0, cor(y_pred, y_test), 0)
    zscore_folds[test_fold] <- atanh(corr_folds[test_fold])*sqrt(length(y_test) - 3) # Fisher transformation
    pval_folds[test_fold] <- ifelse(sd(y_pred) != 0, cor.test(y_pred, y_test)$p.value, runif(1))
  }
  R2_avg <- mean(R2_folds)
  R2_sd <- sd(R2_folds)
  rho_avg <- mean(corr_folds)
  rho_se <- sd(corr_folds)
  rho_avg_squared <- rho_avg**2
  # Stouffer's method for combining z scores.
  zscore_est <- sum(zscore_folds) / sqrt(n_train_test_folds)
  zscore_pval <- 2*pnorm(abs(zscore_est), lower.tail = FALSE)
  # Fisher's method for combining p-values: https://en.wikipedia.org/wiki/Fisher%27s_method
  pval_est <- pchisq(-2 * sum(log(pval_folds)), 2*n_train_test_folds, lower.tail = F)
  list(R2_avg=R2_avg, R2_sd=R2_sd, pval_est=pval_est, rho_avg=rho_avg, rho_se=rho_se, rho_zscore=zscore_est, rho_avg_squared=rho_avg_squared, zscore_pval=zscore_pval)
}

do_covariance <- function(gene_id, cis_gt, rsids, varIDs) {
  model_gt <- cis_gt[,rsids, drop=FALSE]
  colnames(model_gt) <- rsids
  geno_cov <- cov(model_gt)
  geno_cov[lower.tri(geno_cov)] <- NA
  cov_df <- melt(geno_cov, varnames = c("rsid1", "rsid2"), na.rm = TRUE) %>%
              mutate(gene=gene_id) %>%
              select(GENE=gene, RSID1=rsid1, RSID2=rsid2, VALUE=value) %>%
              arrange(GENE, RSID1, RSID2)
  cov_df
}

main <- function(snp_annot_file, gene_annot_file, genotype_file, expression_file,dapg_file,
                 covariates_file, chrom, prefix, maf=0.01, n_folds=10, n_train_test_folds=5,
                 seed=NA, cis_window=1e6, alpha=0.5, null_testing=FALSE, pip_threshold=0.01,filt_cluster=T) { ### added in dapg_file as input to main
  cat("loading gene annotation\n")
  gene_annot <- get_gene_annotation(gene_annot_file, chrom)
  cat("loading gene expression\n") 
  expr_df <- get_gene_expression(expression_file, gene_annot)
  samples <- rownames(expr_df)
  n_samples <- length(samples)
  genes <- colnames(expr_df)
  n_genes <- length(expr_df)
  cat("loading genotype annotation\n") 
  snp_annot <- get_filtered_snp_annot(snp_annot_file)
  cat("loading genotypes\n") 
  gt_df <- get_maf_filtered_genotype(genotype_file, maf, samples,snp_annot)
  covariates_df <- get_covariates(covariates_file, samples)
  PIP_df<-read_dapg(dapg_file) ### get dapg results
  
  # Set seed----
  seed <- ifelse(is.na(seed), sample(1:1000000, 1), seed)
  set.seed(seed)
  registerDoRNG(seed)
  
  # Prepare output data----
  model_summary_file <- prefix %&% '_chr' %&% chrom %&% '_model_summaries.txt'
  model_summary_cols <- c('gene_id', 'gene_name', 'gene_type', 'alpha', 'n_snps_in_window', 'n_snps_in_model', 'lambda_min_mse',
                          'test_R2_avg', 'test_R2_sd', 'cv_R2_avg', 'cv_R2_sd', 'in_sample_R2',
                          'nested_cv_fisher_pval', 'rho_avg', 'rho_se', 'rho_zscore', 'rho_avg_squared', 'zscore_pval',
                          'cv_rho_avg', 'cv_rho_se', 'cv_rho_avg_squared', 'cv_zscore_est', 'cv_zscore_pval', 'cv_pval_est')
  write(model_summary_cols, file = model_summary_file, ncol = 24, sep = '\t')
  
  weights_file <-prefix %&% '_chr' %&% chrom %&% '_weights.txt' ###edited prefix so it now functions like a PLINK prefix
  weights_col <- c('gene_id', 'rsid', 'varID', 'ref', 'alt', 'beta')
  write(weights_col, file = weights_file, ncol = 6, sep = '\t')
  
  tiss_chr_summ_f <- prefix %&% '_chr' %&% chrom %&% '_tiss_chr_summary.txt' ###edited prefix so it now functions like a PLINK prefix
  tiss_chr_summ_col <- c('n_samples', 'chrom', 'cv_seed', 'n_genes')
  tiss_chr_summ <- data.frame(n_samples, chrom, seed, n_genes)
  colnames(tiss_chr_summ) <- tiss_chr_summ_col
  write.table(tiss_chr_summ, file = tiss_chr_summ_f, quote = FALSE, row.names = FALSE, sep = '\t')
  
  covariance_file <-prefix %&% '_chr' %&% chrom %&% '_covariances.txt' ###edited prefix so it now functions like a PLINK prefix
  covariance_col <- c('GENE', 'RSID1', 'RSID2', 'VALUE')
  write(covariance_col, file = covariance_file, ncol = 4, sep = ' ')
  
  # Attempt to build model for each gene----
  for (i in 1:n_genes) {
    ptm <- proc.time()
    cat(i, "/", n_genes, "\n")
    gene <- genes[i]
    gene_name <- gene_annot$gene_name[gene_annot$gene_id == gene]
    gene_type <- get_gene_type(gene_annot, gene)
    coords <- get_gene_coords(gene_annot, gene)
    dapg_gt <- get_dapg_genotype(gt_df, PIP_df,gene, pip_threshold,filt_cluster) ### instead of get_cis_gt now get_dapg_gt
    penalty <- get_dapg_penalty(dapg_gt,PIP_df,gene) ### Also get associated penalty
    if (all(is.na(dapg_gt))) { ### If no snps meet threshold then move on
      # No snps within window for gene.
      model_summary <- c(gene, gene_name, gene_type, alpha, 0, 0, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
      write(model_summary, file = model_summary_file, append = TRUE, ncol = 24, sep = '\t')
      next
    }
    model_summary <- c(gene, gene_name, gene_type, alpha, ncol(dapg_gt), 0, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
    if (ncol(dapg_gt) >= 2) {
      expression_vec <- expr_df[,i]
      # str(expression_vec)
      # print(expression_vec)
      if(!is.na(covariates_df)){ ### check if error handling returned NA, adjust if covariates file exists
        adj_expression <- adjust_for_covariates(expression_vec, covariates_df)
      } else{
        adj_expression <- as.matrix(expression_vec)
      }
      # str(adj_expression)
      if (null_testing)
        adj_expression <- sample(adj_expression)
      perf_measures <- nested_cv_elastic_net_perf(dapg_gt, adj_expression, n_samples, n_train_test_folds, n_folds, alpha,penalty)
      R2_avg <- perf_measures$R2_avg
      R2_sd <- perf_measures$R2_sd
      pval_est <- perf_measures$pval_est
      rho_avg <- perf_measures$rho_avg
      rho_se <- perf_measures$rho_se
      rho_zscore <- perf_measures$rho_zscore
      rho_avg_squared <- perf_measures$rho_avg_squared
      zscore_pval <- perf_measures$zscore_pval
      # Fit on all data
      cv_fold_ids <- generate_fold_ids(length(adj_expression), n_folds)
      #str(dapg_gt)
      #print(class(dapg_gt))
      #str(adj_expression)
      #print(class(adj_expression))
      #str(penalty)
      #print(class(penalty))
      fit <- tryCatch(cv.glmnet(dapg_gt, adj_expression, nfolds = n_folds, alpha = 0.5, type.measure='mse', foldid = cv_fold_ids, keep = TRUE, parallel = TRUE,penalty.factor=penalty),
                      error = function(cond) {message('Error'); message(geterrmessage()); list()}) ### Added in penalty.factor=penalty
      
      if (length(fit) > 0) {
        cv_R2_folds <- rep(0, n_folds)
        cv_corr_folds <- rep(0, n_folds)
        cv_zscore_folds <- rep(0, n_folds)
        cv_pval_folds <- rep(0, n_folds)
        best_lam_ind <- which.min(fit$cvm)
        for (j in 1:n_folds) {
          fold_idxs <- which(cv_fold_ids == j)
          adj_expr_fold_pred <- fit$fit.preval[fold_idxs, best_lam_ind]
          cv_R2_folds[j] <- calc_R2(adj_expression[fold_idxs], adj_expr_fold_pred)
          cv_corr_folds[j] <- ifelse(sd(adj_expr_fold_pred) != 0, cor(adj_expr_fold_pred, adj_expression[fold_idxs]), 0)
          cv_zscore_folds[j] <- atanh(cv_corr_folds[j])*sqrt(length(adj_expression[fold_idxs]) - 3) # Fisher transformation
          cv_pval_folds[j] <- ifelse(sd(adj_expr_fold_pred) != 0, cor.test(adj_expr_fold_pred, adj_expression[fold_idxs])$p.value, runif(1))
        }
        
        cv_R2_avg <- mean(cv_R2_folds)
        cv_R2_sd <- sd(cv_R2_folds)
        adj_expr_pred <- predict(fit, as.matrix(dapg_gt), s = 'lambda.min')
        training_R2 <- calc_R2(adj_expression, adj_expr_pred)
        
        cv_rho_avg <- mean(cv_corr_folds)
        cv_rho_se <- sd(cv_corr_folds)
        cv_rho_avg_squared <- cv_rho_avg**2
        # Stouffer's method for combining z scores.
        cv_zscore_est <- sum(cv_zscore_folds) / sqrt(n_folds)
        cv_zscore_pval <- 2*pnorm(abs(cv_zscore_est), lower.tail = FALSE)
        cv_pval_est <- pchisq(-2 * sum(log(cv_pval_folds)), 2*n_folds, lower.tail = F)
        if (fit$nzero[best_lam_ind] > 0) {
          weights <- fit$glmnet.fit$beta[which(fit$glmnet.fit$beta[,best_lam_ind] != 0), best_lam_ind]
          weighted_snps <- names(fit$glmnet.fit$beta[,best_lam_ind])[which(fit$glmnet.fit$beta[,best_lam_ind] != 0)]
           str(weighted_snps)
          # str(snp_annot)
          weighted_snps_info <- snp_annot %>% filter(rsid %in% weighted_snps) %>% select(rsid, varID, refAllele, effectAllele);###edited column names and which columns are being compared
          ###I edited column names for a bunch of the dplyr filtering here. There should be no major changes to the code base, but my data is not varID so that had to be changed else it would filter out everything.
          weighted_snps_info$gene <- gene;
          str(weighted_snps_info)
          weighted_snps_info <- weighted_snps_info %>%
            merge(data.frame(weights = weights, rsid=weighted_snps), by = 'rsid') %>%
            select(gene, rsid, varID, refAllele, effectAllele, weights)
          str(weighted_snps_info)
          write.table(weighted_snps_info, file = weights_file, append = TRUE, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = '\t')
          covariance_df <- do_covariance(gene, dapg_gt, weighted_snps_info$rsid, weighted_snps_info$varID)
          write.table(covariance_df, file = covariance_file, append = TRUE, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = " ")
          model_summary <- c(gene, gene_name, gene_type, alpha, ncol(dapg_gt), fit$nzero[best_lam_ind], fit$lambda[best_lam_ind], R2_avg, R2_sd, cv_R2_avg, cv_R2_sd, training_R2, pval_est,
                             rho_avg, rho_se, rho_zscore, rho_avg_squared, zscore_pval, cv_rho_avg, cv_rho_se, cv_rho_avg_squared, cv_zscore_est, cv_zscore_pval, cv_pval_est)
        } else {
          model_summary <- c(gene, gene_name, gene_type, alpha, ncol(dapg_gt), 0, fit$lambda[best_lam_ind], R2_avg, R2_sd,
                             cv_R2_avg, cv_R2_sd, training_R2, pval_est, rho_avg, rho_se, rho_zscore, rho_avg_squared, zscore_pval,
                             cv_rho_avg, cv_rho_se, cv_rho_avg_squared, cv_zscore_est, cv_zscore_pval, cv_pval_est)
        }
      } else {
        model_summary <- c(gene, gene_name, gene_type, alpha, ncol(dapg_gt), 0, NA, R2_avg, R2_sd, NA, NA, NA, pval_est, rho_avg, rho_se, rho_zscore, rho_avg_squared, zscore_pval,
                           NA, NA, NA, NA, NA, NA)
      }
    }
    write(model_summary, file = model_summary_file, append = TRUE, ncol = 24, sep = '\t')
    proc.time() - ptm
  }
}### We should consider adding in having it write the penalty factor of each snp to the model summary
