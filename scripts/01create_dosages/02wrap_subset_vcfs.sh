#!/bin/bash
#PBS -S /bin/bash
#PBS -l nodes=1:ppn=1
#PBS -o /home/ryan/logs/${PBS_JOBNAME}.o${PBS_JOBID}.log
#PBS -e /home/ryan/logs/${PBS_JOBNAME}.e${PBS_JOBID}.err
#PBS -l walltime=150:00:00
#PBS -l mem=16gb

/usr/local/bin/bcftools view -S /home/ryan/topmed/proteome/sample_lists/${pop}_sidno_id_list.txt \
-O z \
-o /home/ryan/topmed/proteome/${pop}/genotypes/vcfs/proteome_${pop}_chr${chr}.vcf.gz \
/home/wheelerlab3/Data/TOPMed/MESA_TOPMED_Imputation/${pop}/chr${chr}.dose.vcf.gz
