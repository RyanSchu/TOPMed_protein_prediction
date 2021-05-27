#!/bin/bash

#$ -S /bin/bash
#$ -pe mpi 1

echo $pop

/home/rschubert1/anaconda3/bin/python3.8 /scratch/rschubert1/MetaXcan/software/PrediXcan.py \
--throw \
--text_genotypes /data/rschubert1/topmed_proteome_ashg_prep/liftover/noheader_hg38chr*.maf0.01.R20.8.dosage.txt.gz \
--text_sample_ids /data/rschubert1/topmed_proteome_ashg_prep/liftover/prediXcan_samples.txt \
--model_db_path /scratch/rschubert1/TOPMed_Proteome/06Elastic_net/PAV_modeling/db_out/${pop}_PCAIR_baseline_PAV_adj_models_rho0.1_zpval0.05.db \
--prediction_output /home/rschubert1/scratch/TOPMed_Proteome/10INTERVAL_prediction/PAV_adj_predictions/${pop}_PCAIR_baseline_PAV_adj_models_rho0.1_zpval0.05_prediction.txt \
--prediction_summary_output /home/rschubert1/scratch/TOPMed_Proteome/10INTERVAL_prediction/PAV_adj_predictions/${pop}_PCAIR_baseline_PAV_adj_models_rho0.1_zpval0.05_summary.txt


#echo /home/rschubert1/data/topmed_proteome_ashg_prep/imputed_INTERVAL/${pop}_baseline_unfiltered_model_INTERVAL_predicted_expression.txt
