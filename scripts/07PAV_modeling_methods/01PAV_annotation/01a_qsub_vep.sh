#!/bin/bash
#PBS -S /bin/bash
#PBS -l walltime=36:00:00
#PBS -l nodes=1:ppn=2
#PBS -l mem=16gb
#PBS -e /home/ryan/logs/${PBS_JOBNAME}.err
#PBS -o /home/ryan/logs/${PBS_JOBNAME}.out

chr=22
pop=CAU
export PERL5LIB=/home/ryan/anaconda3/bin/perl

#perl -V

/home/ryan/software/ensembl-vep/vep \
-i /home/wheelerlab3/Data/TOPMed/MESA_TOPMED_Imputation/${pop}/chr${chr}.dose.vcf.gz \
--cache \
--per_gene \
--assembly GRCh37 \
-o /home/ryan/software/ensembl-vep/${pop}_anno/test_chr${chr}.dose.anno.txt \
--force_overwrite \
--show_cache_info
