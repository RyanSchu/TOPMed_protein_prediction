#!/bin/bash

#$ -S /bin/bash
#$ -pe mpi 1

echo $db

/home/rschubert1/anaconda3/bin/python3.8 /scratch/rschubert1/MetaXcan/software/PrediXcan.py \
--throw \
--text_genotypes /data/rschubert1/TOPMED_Proteome/INTERVAL/liftover/hg38chr*.maf0.01.R20.8.dosage_TRUE_SUBSET.txt.gz \
--text_sample_ids /home/rschubert1/data/TOPMED_Proteome/INTERVAL/liftover/prediXcan_samples.txt \
--model_db_path ${db} \
--prediction_output /home/rschubert1/scratch/TOPMed_Proteome/10INTERVAL_prediction/true_subset_predictions/${db:46:(-3)}_predicted_expression.txt \
--prediction_summary_output /home/rschubert1/scratch/TOPMed_Proteome/10INTERVAL_prediction/true_subset_predictions/${db:46:(-3)}_summary.txt

#echo /home/rschubert1/data/topmed_proteome_ashg_prep/imputed_INTERVAL/${pop}_baseline_unfiltered_model_INTERVAL_predicted_expression.txt
