#!/bin/bash
pop=$1
Rscript /home/ryan/TOPMed_Proteome/03adjust_expression/03reformat_protein_ids.R \
--proteins /home/ryan/topmed/proteome/${pop}/PCAIR_modeling/03adjusted_expression/Proteome_TOPMed_${pop}_ln_adjAgeSex_mean_rank-inverse_adj10PCs_MEQTL.txt.gz \
--annotation /home/ryan/topmed/proteome/dapg/Aptamer_annotation.txt.gz \
--out /home/ryan/topmed/proteome/${pop}/PCAIR_modeling/03adjusted_expression/Proteome_TOPMed_${pop}_ln_adjAgeSex_mean_rank-inverse_adj10PCair_MEQTL_reformatted_ids.txt
gzip /home/ryan/topmed/proteome/${pop}/PCAIR_modeling/03adjusted_expression/Proteome_TOPMed_${pop}_ln_adjAgeSex_mean_rank-inverse_adj10PCair_MEQTL_reformatted_ids.txt
