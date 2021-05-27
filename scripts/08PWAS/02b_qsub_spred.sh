#!/bin/bash

for p in CHN AFA CAU ALL HIS
do
	for i in /home/rschubert1/data/ukbb_sum_stats/spred_input/Whradjbmi.giant-ukbb.meta-analysis.*males.23May2018-.HapMap2_only.hg38.txt.gz
	do
		for m in /home/rschubert1/data/TOPMED_Proteome/${p}/PCAIR_modeling/06Elastic_net/dbs_out/*rho0.1_zpval0.05.db
		do
			qsub -v pop=${p},file=${i},model=${m} -N SPred_${p} /data/rschubert1/ukbb_sum_stats/scripts/SPrediXcan_topmed.sh
		done
	done
done
