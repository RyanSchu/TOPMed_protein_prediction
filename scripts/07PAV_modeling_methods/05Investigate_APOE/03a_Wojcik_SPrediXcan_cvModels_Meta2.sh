cd /home/egeoffroy/MetaXcan_v3/software/
gwas_SS=$1
name=$2
population=('AFA' 'ALL')

for pop in ${population[@]};
do
	/home/rschubert1/anaconda3/bin/python3.8 ./SPrediXcan.py \
	--model_db_path /scratch/rschubert1/TOPMed_Proteome/investigate_APOE/model_PAV_adj_APOE/db/${pop}_PCAIR_baseline_PAV_adj_models_APOE_unfiltered.db \
	--covariance /scratch/rschubert1/TOPMed_Proteome/investigate_APOE/model_PAV_adj_APOE/db/${pop}_PCAIR_baseline_PAV_adj_models_APOE_unfiltered_covariances.txt \
	--keep_non_rsid \
	--gwas_file ${gwas_SS}  \
	--snp_column SNP_hg38 \
	--effect_allele_column Effect-allele \
	--non_effect_allele_column Other-allele \
	--beta_column Beta \
	--pvalue_column P-val \
	--output_file /home/rschubert1/scratch/TOPMed_Proteome/investigate_APOE/PAGE_spred/${name}_${pop}_APOE_PAV_adj.csv
done
