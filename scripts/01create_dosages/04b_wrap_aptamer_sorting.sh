#!/bin/bash

pop=$1

Rscript /home/ryan/Useful_things/utility_scripts/match_order.R \
--data /home/wheelerlab3/topmed/adj_omics_traits/MatrixEQTL_format/Proteome_TOPMed_${pop}_ln_adjAgeSex_mean_rank-inverse_MEQTL.txt \
--match /home/ryan/topmed/proteome/${pop}/genotypes/dosages/chr22.maf0.01.R20.8.dosage.txt.gz \
--skip.col.data 1 \
--skip.col.match 6 \
--out /home/ryan/topmed/proteome/${pop}/pQTL/sorted_input/Proteome_TOPMed_${pop}_ln_adjAgeSex_mean_rank-inverse_MEQTL_sorted.txt
