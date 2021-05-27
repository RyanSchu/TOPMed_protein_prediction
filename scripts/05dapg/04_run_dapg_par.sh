#!/bin/bash

#PBS -S /bin/bash
#PBS -l nodes=1:ppn=40
#PBS -o /home/ryan/logs/${PBS_JOBNAME}.o${PBS_JOBID}.log
#PBS -e /home/ryan/logs/${PBS_JOBNAME}.e${PBS_JOBID}.err
#PBS -l walltime=150:00:00
#PBS -l mem=100gb

#pop=CAU
#tiss=Tcell
proteome_anno=/home/ryan/topmed/proteome/dapg/Aptamer_annotation.txt.gz
#dat_dir=/home/ryan/topmed/expression/${pop}/${tiss}/dapg_in/dat
#priordir=/home/ryan/topmed/expression/CAU/Tcell/dapg_in/redo_priors
#outputdir=/home/ryan/topmed/expression/${pop}/${tiss}/dapg_out/

if [ ! -d  ${outputdir} ]
then
	mkdir ${outputdir}
fi

run_dapg() {
	pop=$1
	arr=($2)
	line=${arr[0]}
	dat_dir=/home/ryan/topmed/proteome/${pop}/PCAIR_modeling/05dapg/dat_files
	priordir=/home/ryan/topmed/proteome/${pop}/PCAIR_modeling/05dapg/torus_out_priors
	outputdir=/home/ryan/topmed/proteome/${pop}/PCAIR_modeling/05dapg/dapg_out
	if [ -e "${dat_dir}/${line}.dat" ]
        then
                echo "${dat_dir}/${line}.dat"
                echo "${priordir}/${line}.prior"
                /usr/local/bin/dap-g -d ${dat_dir}/${line}.dat -p ${priordir}/${line}.prior -t 1  -o ${outputdir}/${line}.fm.out -l ${outputdir}/${line}.log
	fi
}


export -f run_dapg

for i in CHN AFA HIS CAU ALL
do
	zcat $proteome_anno | parallel -j 30 run_dapg ${i} {}
done
