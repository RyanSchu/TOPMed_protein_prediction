## make PC-Air PCs

args = commandArgs(trailingOnly=TRUE)
"%&%" = function(a,b) paste(a,b,sep="")

library(GENESIS)
library(SNPRelate)


snpgdsBED2GDS(bed.fn="/home/ryan/HAPMAP3_hg19/hg38/topmed_id_format/hapmap.hg38.alternate_topmed_id.bed",
              bim.fn="/home/ryan/HAPMAP3_hg19/hg38/topmed_id_format/hapmap.hg38.alternate_topmed_id.bim",
              fam.fn="/home/ryan/HAPMAP3_hg19/hg38/topmed_id_format/hapmap.hg38.alternate_topmed_id.fam",
              out.gdsfn="/home/ryan/HAPMAP3_hg19/hg38/topmed_id_format/hapmap.hg38.alternate_topmed_id.gds")

