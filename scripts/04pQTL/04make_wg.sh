#!/bin/bash

for i in AFA CHN HIS ALL CAU
do
	cat /home/ryan/topmed/proteome/${i}/PCAIR_modeling/04pQTL/cis_eQTLs_${i}_chr1_all_cis.txt > /home/ryan/topmed/proteome/${i}/PCAIR_modeling/04pQTL/cis_eQTLs_${i}_WG_all_cis.txt
	gzip /home/ryan/topmed/proteome/${i}/PCAIR_modeling/04pQTL/cis_eQTLs_${i}_chr1_all_cis.txt
	for j in {2..22}
	do
		echo $i $j
		tail -n +2 /home/ryan/topmed/proteome/${i}/PCAIR_modeling/04pQTL/cis_eQTLs_${i}_chr${j}_all_cis.txt >> /home/ryan/topmed/proteome/${i}/PCAIR_modeling/04pQTL/cis_eQTLs_${i}_WG_all_cis.txt
		gzip /home/ryan/topmed/proteome/${i}/PCAIR_modeling/04pQTL/cis_eQTLs_${i}_chr${j}_all_cis.txt
	done
	gzip /home/ryan/topmed/proteome/${i}/PCAIR_modeling/04pQTL/cis_eQTLs_${i}_WG_all_cis.txt
done
