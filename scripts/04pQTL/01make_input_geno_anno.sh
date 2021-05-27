#!/bin/bash
#PBS -S /bin/bash
#PBS -l nodes=1:ppn=1
#PBS -o /home/ryan/logs/${PBS_JOBNAME}.o${PBS_JOBID}.log
#PBS -e /home/ryan/logs/${PBS_JOBNAME}.e${PBS_JOBID}.err
#PBS -l walltime=150:00:00
#PBS -l mem=16gb


pop=ALL
for i in {22..22}
do
	echo $pop $i
	dosage=/home/ryan/topmed/proteome/${pop}/PCAIR_modeling/01genotypes/dosages/hg38chr${i}.maf0.01.R20.8.dosage.txt.gz
	outgeno=/home/ryan/topmed/proteome/${pop}/PCAIR_modeling/01genotypes/preddb_input/uniq_pred_db_hg38.chr${i}.maf0.01.R20.8.geno.txt
	outanno=/home/ryan/topmed/proteome/${pop}/PCAIR_modeling/01genotypes/preddb_input//uniq_pred_db_hg38.chr${i}.maf0.01.R20.8.anno.txt
	MEQTLanno=/home/ryan/topmed/proteome/${pop}/PCAIR_modeling/01genotypes/MEQTL_input/MEQTL_annotation_chr${i}.txt
	tmp=/home/ryan/topmed/proteome/${pop}/genotypes/preddb_input_PAV_filtered/tmp_degenerate_list_chr${i}.txt
	#make the genotype files
	zcat ${dosage} | cut -f 2 -d " " | sort | uniq -d > ${tmp}
	zcat ${dosage} | cut -f 2,7- -d " " | grep -v -F -f ${tmp} > ${outgeno}
	#make the annot file
	echo "chr pos varID refAllele effectAllele rsid" > ${outanno}
	zcat ${dosage} | tail -n +2 | grep -v -F -f ${tmp} | awk '{gsub("chr","",$1);print $1 FS $3 FS $1 "_" $3 "_" $4 "_" $5 "_b38" FS $4 FS $5 FS $2}' >> ${outanno}
	awk '{print $6 FS $1 FS $2}' ${outanno} > ${MEQTLanno}
	gzip -f ${outanno}
	gzip -f ${outgeno}
	gzip -f ${MEQTLanno}
	rm ${tmp}
done
