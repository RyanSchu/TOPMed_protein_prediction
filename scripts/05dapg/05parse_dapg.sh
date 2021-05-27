#!/bin/bash


dapg_dir=$1
output_dir=$2
tagDefault=summary
tag=$3
tag=${tag:=$tagDefault}

if [ -e ${output_dir}/${tag}_snps.txt ]
then
	rm ${output_dir}/${tag}_snps.txt
fi

if [ -e ${output_dir}/${tag}_clusters.txt ]
then
	rm ${output_dir}/${tag}_clusters.txt
fi

for i in ${dapg_dir}/*.fm.out
do
	grep -F "((" ${i} | awk -v txt=$( basename ${i} ) '{if ($5 != -1) print $0 FS txt}' >> ${output_dir}/${tag}_snps.txt
	grep -F "((" ${i} | awk -v txt=$( basename ${i} ) '{print $0 FS txt}' >> ${output_dir}/${tag}_snps_include_unassigned.txt
	grep -F "{" ${i} | awk -v txt=$( basename ${i} ) '{print $1 FS $2 FS $3 FS $4 FS txt}' >> ${output_dir}/${tag}_clusters.txt
done
