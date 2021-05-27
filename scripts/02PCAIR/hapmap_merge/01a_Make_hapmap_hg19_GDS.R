## make PC-Air PCs

args = commandArgs(trailingOnly=TRUE)
"%&%" = function(a,b) paste(a,b,sep="")

library(GENESIS)
library(SNPRelate)


snpgdsBED2GDS(bed.fn="/home/ryan/HAPMAP3_hg19/topmed_id_format/HM3_ASN_CEU_YRI_Unrelated_hg19_noAmbig.bed",
              bim.fn="/home/ryan/HAPMAP3_hg19/topmed_id_format/HM3_ASN_CEU_YRI_Unrelated_hg19_noAmbig.bim",
              fam.fn="/home/ryan/HAPMAP3_hg19/topmed_id_format/HM3_ASN_CEU_YRI_Unrelated_hg19_noAmbig.fam",
              out.gdsfn="/home/ryan/HAPMAP3_hg19/topmed_id_format/HM3_ASN_CEU_YRI_Unrelated_hg19_noAmbig.gds")

