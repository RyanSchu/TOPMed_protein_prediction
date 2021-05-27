#!/bin/bash
#PBS -S /bin/bash
#PBS -l walltime=36:00:00
#PBS -l nodes=1:ppn=2
#PBS -l mem=16gb
#PBS -e /home/ryan/logs/${PBS_JOBNAME}.err
#PBS -o /home/ryan/logs/${PBS_JOBNAME}.out
#pop=$1
#chr=$2

if [ $pop == "CAU" ]
then
	input=/home/ryan/software/ensembl-vep/${pop}_anno/hg19_chr${chr}.dose.anno.txt
else
	input=/home/ryan/software/ensembl-vep/${pop}_anno/chr${chr}.dose.anno.txt
fi

/home/ryan/software/ensembl-vep/filter_vep \
-i ${input} \
--filter "Consequence in coding_sequence_variant,frameshift_variant,inframe_deletion,inframe_insertion,missense_variant,protein_altering_variant,splice_acceptor_variant,splice_donor_variant,splice_region_variant,start_lost,stop_gained,stop_lost" \
-o /home/ryan/software/ensembl-vep/${pop}_anno/protein_altering_snps_chr${chr}.txt
