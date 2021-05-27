#!/bin/bash

for po in ALL
do
	for chr in {1..22}
	do
		for pi in 0 0.001 0.01
		do
			for cls in T F
			do
				qsub -v pop=${po},chrom=${chr},pip=${pi},clus=${cls} -N ${po}_${chr}_${pi}_${cls}_PAV_filtered 02f_TOPMED_dapg_qsub_PAV_filtered.pbs
			done
		done
	done
done
