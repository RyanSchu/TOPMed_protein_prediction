#!/bin/bash

#$ -S /bin/bash
#$ -pe mpi 1

echo $pop $file

/home/rschubert1/anaconda3/bin/python3.8 /scratch/rschubert1/MetaXcan/software/SPrediXcan.py \
--model_db_path ${model} \
--covariance ${model::-3}_covariances.txt \
--gwas_file ${file} \
--keep_non_rsid \
--snp_column SNP_hg38 \
--effect_allele_column alt \
--non_effect_allele_column ref \
--beta_column beta \
--pvalue_column pval \
--output_file ${file::-7}_Spred_${model:79:(-3)}

mv ${file::-7}_Spred_${model:79:(-3)} /home/rschubert1/data/ukbb_sum_stats/spred_out/
