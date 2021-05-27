#!/bin/bash
pop=$1

Rscript /home/ryan/software/Useful_things/utility_scripts/match_order.R \
--data /home/ryan/topmed/proteome/${pop}/PCAIR_modeling/03adjusted_expression/Proteome_TOPMed_${pop}_ln_adjAgeSex_mean_rank-inverse_adj10PCair_MEQTL_reformatted_ids.txt.gz \
--match /home/ryan/topmed/proteome/${pop}/PCAIR_modeling/01genotypes/preddb_input/uniq_pred_db_hg38.chr22.maf0.01.R20.8.geno.txt.gz \
--skip.col.data 1 \
--skip.col.match 1 \
--out /home/ryan/topmed/proteome/${pop}/PCAIR_modeling/03adjusted_expression/Proteome_TOPMed_${pop}_ln_adjAgeSex_mean_rank-inverse_adj10PCair_MEQTL_reformatted_ids_sorted.txt
gzip -f /home/ryan/topmed/proteome/${pop}/PCAIR_modeling/03adjusted_expression/Proteome_TOPMed_${pop}_ln_adjAgeSex_mean_rank-inverse_adj10PCair_MEQTL_reformatted_ids_sorted.txt
