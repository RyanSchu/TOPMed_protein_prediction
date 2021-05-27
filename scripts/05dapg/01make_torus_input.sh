#!/bin/bash
pop=$1
echo $pop
zcat /home/ryan/topmed/proteome/${pop}/PCAIR_modeling/04pQTL/cis_eQTLs_${pop}_WG_all_cis.txt.gz | awk '{print $1 "\t" $2 "\t" $6 "\t" $3 "\t" $4}' > /home/ryan/topmed/proteome/${pop}/PCAIR_modeling/05dapg/cis_eQTLs_${pop}_all_cis_unpruned_TORUS_input.txt
tail -n +2 /home/ryan/topmed/proteome/${pop}/PCAIR_modeling/05dapg/cis_eQTLs_${pop}_all_cis_unpruned_TORUS_input.txt | awk '{split($1,a,":");gsub("chr","",a[1]);print $1 "\t" a[1] "\t" a[2]}' | sort -gk2,2 -gk3,3 | uniq > /home/ryan/topmed/proteome/${pop}/PCAIR_modeling/05dapg/${pop}_snp_annotation.txt
gzip -f /home/ryan/topmed/proteome/${pop}/PCAIR_modeling/05dapg/${pop}_snp_annotation.txt
gzip -f /home/ryan/topmed/proteome/${pop}/PCAIR_modeling/05dapg/cis_eQTLs_${pop}_all_cis_unpruned_TORUS_input.txt
