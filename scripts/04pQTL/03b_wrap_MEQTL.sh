#!/bin/bash
#PBS -S /bin/bash
#PBS -l nodes=1:ppn=1
#PBS -o /home/ryan/logs/${PBS_JOBNAME}.o${PBS_JOBID}.log
#PBS -e /home/ryan/logs/${PBS_JOBNAME}.e${PBS_JOBID}.err
#PBS -l walltime=150:00:00
#PBS -l mem=100gb

pop=$1
chr=$2

if [ $pop == "ALL" ]
then
	expression=/home/ryan/topmed/proteome/ALL/PCAIR_modeling/03adjusted_expression/Proteome_TOPMed_ALL_ln_adjAgeSex_mean_rank-inverse_adj10PCair_MEQTL_reformatted_ids_sorted.txt.gz
else
	expression=/home/ryan/topmed/proteome/${pop}/PCAIR_modeling/03adjusted_expression/Proteome_TOPMed_${pop}_ln_adjAgeSex_mean_rank-inverse_adj10PCair_MEQTL_reformatted_ids.txt.gz
fi

Rscript /home/ryan/TOPMed_Proteome/04pQTL/03a_MatrixEQTL.R \
--snpgenotype /home/ryan/topmed/proteome/${pop}/PCAIR_modeling/01genotypes/preddb_input/uniq_pred_db_hg38.chr${chr}.maf0.01.R20.8.geno.txt.gz \
--snplocation /home/ryan/topmed/proteome/${pop}/PCAIR_modeling/01genotypes/MEQTL_input/MEQTL_annotation_chr${chr}.txt.gz \
--geneexpression ${expression} \
--genelocation /home/ryan/topmed/proteome/ALL/PCAIR_ALL/02pQTL/aptamer_annotation_v2.txt.gz \
--tag ${pop}_chr${chr}_all_cis \
--outputdir /home/ryan/topmed/proteome/${pop}/PCAIR_modeling/04pQTL/ \
--cis 1 \
--trans 0 \
--window 1000000
