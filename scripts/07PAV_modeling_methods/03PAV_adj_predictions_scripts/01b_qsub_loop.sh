#!/bin/bash

for p in ALL AFA CHN HIS CAU
do
	qsub -v pop=${p} -N ${p}_PAV_adj_prediction 01impute_db.sh
done
