
for i in AFA CHN CAU HIS
do
	zcat /data/rschubert1/TOPMED_Proteome/${i}/PCAIR_modeling/03adjusted_expression/Proteome_TOPMed_${i}_ln_adjAgeSex_mean_rank-inverse_adj10PCair_MEQTL_reformatted_ids.txt.gz | head -n 1 > /home/rschubert1/scratch/TOPMed_Proteome/investigate_APOE/${i}_expression.txt
	zgrep -f /scratch/rschubert1/TOPMed_Proteome/investigate_APOE/APOE_aptamers.txt /data/rschubert1/TOPMED_Proteome/${i}/PCAIR_modeling/03adjusted_expression/Proteome_TOPMed_${i}_ln_adjAgeSex_mean_rank-inverse_adj10PCair_MEQTL_reformatted_ids.txt.gz >> /home/rschubert1/scratch/TOPMed_Proteome/investigate_APOE/${i}_expression.txt
done

zcat /data/rschubert1/TOPMED_Proteome/ALL/PCAIR_modeling/03adjusted_expression/Proteome_TOPMed_ALL_ln_adjAgeSex_mean_rank-inverse_adj10PCair_MEQTL_reformatted_ids_sorted.txt.gz | head -n 1 > /home/rschubert1/scratch/TOPMed_Proteome/investigate_APOE/ALL_expression.txt
zgrep -f /scratch/rschubert1/TOPMed_Proteome/investigate_APOE/APOE_aptamers.txt /data/rschubert1/TOPMED_Proteome/ALL/PCAIR_modeling/03adjusted_expression/Proteome_TOPMed_ALL_ln_adjAgeSex_mean_rank-inverse_adj10PCair_MEQTL_reformatted_ids_sorted.txt.gz >> /home/rschubert1/scratch/TOPMed_Proteome/investigate_APOE/ALL_expression.txt
