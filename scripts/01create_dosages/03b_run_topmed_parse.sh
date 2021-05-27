#!/bin/bash
#PBS -S /bin/bash
#PBS -l nodes=1:ppn=1
#PBS -o /home/ryan/logs/${PBS_JOBNAME}.o${PBS_JOBID}.log
#PBS -e /home/ryan/logs/${PBS_JOBNAME}.e${PBS_JOBID}.err
#PBS -l walltime=150:00:00
#PBS -l mem=16gb

pop=CAU
python2 /home/ryan/Imputation/topmed.py \
-i /home/ryan/topmed/proteome/${pop}/genotypes/vcfs/proteome_${pop}_ \
-c ${chr} \
--cpos \
--outdir /home/ryan/topmed/proteome/${pop}/genotypes/dos_from_lauren
