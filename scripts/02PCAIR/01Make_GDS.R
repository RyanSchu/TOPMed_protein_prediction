## make PC-Air PCs

args = commandArgs(trailingOnly=TRUE)
"%&%" = function(a,b) paste(a,b,sep="")

library(GENESIS)
library(SNPRelate)

if (args[1] == "CAU") {
	dir<-"QC_lauren"
} else {
	dir<-"QC"
}

snpgdsBED2GDS(bed.fn="/home/ryan/topmed/proteome/" %&% args[1] %&%  "/genotypes/" %&% dir %&%  "/missingness_hwe_steps/00autosome.bed",
              bim.fn="/home/ryan/topmed/proteome/" %&% args[1] %&%  "/genotypes/" %&% dir %&%  "/missingness_hwe_steps/00autosome.bim",
              fam.fn="/home/ryan/topmed/proteome/" %&% args[1] %&%  "/genotypes/" %&% dir %&%  "/missingness_hwe_steps/00autosome.fam",
              out.gdsfn="/home/ryan/topmed/proteome/" %&% args[1] %&%  "/PCAIR/00autosome.gds")
