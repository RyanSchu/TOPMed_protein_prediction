#!/bin/bash
#	mkdir /home/ryan/topmed/proteome/${i}/pQTL/dapg/redo_${i}_priors
pop=$1
torus -d /home/ryan/topmed/proteome/${pop}/PCAIR_modeling/05dapg/cis_eQTLs_${pop}_all_cis_unpruned_TORUS_input.txt.gz \
-smap /home/ryan/topmed/proteome/${pop}/PCAIR_modeling/05dapg/${pop}_snp_annotation.txt.gz \
-gmap /home/ryan/topmed/proteome/dapg/Aptamer_annotation.txt.gz \
-dump_prior /home/ryan/topmed/proteome/${pop}/PCAIR_modeling/05dapg/torus_out_priors
