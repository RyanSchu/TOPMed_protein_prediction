#!/bin/bash
for pop in ALL
do
	Rscript /home/ryan/TOPMed_Proteome/04pQTL/07plot_cis_trans_manhattan.R \
	--cis /home/ryan/topmed/proteome/${pop}/PCAIR_modeling/04pQTL/cis_eQTLs_${pop}_WG_all_cis.txt.gz \
	--trans /home/ryan/topmed/proteome/${pop}/PCAIR_modeling/04pQTL/04b_trans_QTLs/${pop}true_trans_pQTL_1e-4.txt.gz \
	--title ${pop} \
	--out /home/ryan/topmed/proteome/${pop}/PCAIR_modeling/04pQTL/${pop}_manhattan_1e-4
done
