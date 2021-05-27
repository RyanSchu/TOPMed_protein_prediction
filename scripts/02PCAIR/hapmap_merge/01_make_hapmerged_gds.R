## make PC-Air PCs

args = commandArgs(trailingOnly=TRUE)
"%&%" = function(a,b) paste(a,b,sep="")

library(GENESIS)
library(SNPRelate)

if (args[1] == "CAU"){
  dir<-"QC_lauren"
} else {
  dir<-"QC"
}

snpgdsBED2GDS(bed.fn="/home/ryan/topmed/proteome/" %&% args[1] %&% "/genotypes/" %&% dir %&% "/hapmap/07final_LDpruned.bed",
              bim.fn="/home/ryan/topmed/proteome/" %&% args[1] %&% "/genotypes/" %&% dir %&% "/hapmap/07final_LDpruned.bim",
              fam.fn="/home/ryan/topmed/proteome/" %&% args[1] %&% "/genotypes/" %&% dir %&% "/hapmap/07final_LDpruned.fam",
              out.gdsfn="/home/ryan/topmed/proteome/" %&% args[1] %&% "/genotypes/" %&% dir %&% "/hapmap/07final_LDpruned.gds")

