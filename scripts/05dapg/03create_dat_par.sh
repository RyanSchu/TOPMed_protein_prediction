#!/bin/bash
#PBS -S /bin/bash
#PBS -l nodes=1:ppn=40
#PBS -o /home/ryan/logs/${PBS_JOBNAME}.o${PBS_JOBID}.log
#PBS -e /home/ryan/logs/${PBS_JOBNAME}.e${PBS_JOBID}.err
#PBS -l walltime=150:00:00
#PBS -l mem=100gb

proteome_anno=/home/ryan/topmed/proteome/dapg/Aptamer_annotation.txt.gz

process_files() {
	pop=$1
	line=($2)
	name=${line[0]}
	chr=${line[1]}

	if [ $pop == "ALL" ]
	then
		proteome_file=/home/ryan/topmed/proteome/${pop}/PCAIR_modeling/03adjusted_expression/Proteome_TOPMed_${pop}_ln_adjAgeSex_mean_rank-inverse_adj10PCair_MEQTL_reformatted_ids_sorted.txt.gz
	else
		proteome_file=/home/ryan/topmed/proteome/${pop}/PCAIR_modeling/03adjusted_expression/Proteome_TOPMed_${pop}_ln_adjAgeSex_mean_rank-inverse_adj10PCair_MEQTL_reformatted_ids.txt.gz
	fi

	genotype_file=/home/ryan/topmed/proteome/${pop}/PCAIR_modeling/01genotypes/preddb_input/uniq_pred_db_hg38.chr${chr}.maf0.01.R20.8.geno.txt.gz

	outdir=/home/ryan/topmed/proteome/${pop}/PCAIR_modeling/05dapg/dat_files
	prior_dir=/home/ryan/topmed/proteome/${pop}/PCAIR_modeling/05dapg/torus_out_priors
	if [ ! -d ${outdir} ]
	then
		mkdir ${outdir}
	fi
	echo "searching $name"
	if [ -e ${prior_dir}/${name}.prior ]
	then
		echo "processing $name"
		zgrep -F ${name} ${proteome_file} | awk -v group=${pop} '{$1="pheno" FS $1 FS group; print }' > ${outdir}/${name}.dat ;
		zgrep -F -f <( cut -f 1 -d ' ' ${prior_dir}/${name}.prior ) ${genotype_file} | sed "s/,/ /g" |awk -v group=${pop} '{$1="geno" FS $1 FS group ;print }' >> ${outdir}/${name}.dat;
	fi
}

export -f process_files

for i in CHN AFA HIS CAU ALL
do
	echo "\n ${i} \n"
	zcat $proteome_anno | parallel -j 30 process_files ${i} {}
done
