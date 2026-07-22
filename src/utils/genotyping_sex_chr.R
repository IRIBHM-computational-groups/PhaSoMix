library(data.table)
library(dplyr)

source("src/utils/writevcf.r")

# --------------------------------------------------------------
# Paths
# --------------------------------------------------------------
path_vcf_MT <- "rawdata/1kGp/rawdata/ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf.gz"
path_vcf_Y  <- "rawdata/1kGp/rawdata/ALL.chrY.phase3_integrated_v2b.20130502.genotypes.vcf.gz"
path_vcf_X  <- "rawdata/1kGp/rawdata/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1c.20130502.genotypes.vcf.gz"

info_path         <- "output/local_ancestry/prediction/high_coverage/IS_corrected/input_path_F1_like_mono_x_admixed_gnomix_PCAWG_POP_no_AMR.tsv"
path_allele_count <- "rawdata/PCAWG/allelecount/"
output_vcf_PCAWG      <- "rawdata/PCAWG/"
output_vcf_1kGp_PCAWG <- "rawdata/1kGp_PCAWG/rawdata/"

vcf_headers <- list(
   MT = "rawdata/vcf_headers/1kGp_hg19_header_mtDNA.txt",
   Y  = "rawdata/vcf_headers/1kGp_hg19_header_chrY.txt",
   X  = "rawdata/vcf_headers/1kGp_hg19_header_chrX.txt"
)

# --------------------------------------------------------------
# Functions
# --------------------------------------------------------------
Getgenotyping <- function(vcf, allel) {
   colnames(allel)[1:2] <- colnames(vcf)[1:2]
   vcf <- merge(vcf[, 1:9], allel)
   for (base in c("A", "C", "G", "T")) {
      vcf$read_ALT[vcf$ALT == base] <- vcf[[paste0("Count_", base)]][vcf$ALT == base]
      vcf$read_REF[vcf$REF == base] <- vcf[[paste0("Count_", base)]][vcf$REF == base]
   }
   vcf$read_ALT   <- as.numeric(vcf$read_ALT)
   vcf$read_REF   <- as.numeric(vcf$read_REF)
   vcf$Good_depth <- as.numeric(vcf$Good_depth)
   vcf$genotype   <- ifelse(vcf$read_ALT > vcf$read_REF & vcf$Good_depth >= 10, 1,
                     ifelse(vcf$read_REF > vcf$read_ALT & vcf$Good_depth >= 10, 0, "."))
   vcf[, c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","genotype")]
}

build_vcf <- function(vcf_ref, samples, chr_tag) {
   vcf_list <- lapply(samples, function(pid) {
      ac_file <- file.path(path_allele_count, paste0(pid, "_alleCount_chr", chr_tag, ".txt"))
      if (!file.exists(ac_file)) return(NULL)
      alFreq  <- fread(ac_file, data.table = FALSE)
      new_vcf <- Getgenotyping(vcf_ref, alFreq)
      colnames(new_vcf)[ncol(new_vcf)] <- pid
      new_vcf
   })
   vcf_list <- Filter(Negate(is.null), vcf_list)
   if (length(vcf_list) == 0) return(NULL)
   Reduce(function(x, y) cbind(x, y[, ncol(y), drop = FALSE]), vcf_list)
}

write_both_vcf <- function(vcf_pcawg, vcf_ref, tag, chr_tag, suffix) {
   vcf_1kgp <- merge(vcf_ref, vcf_pcawg)
   writevcf(vcf_pcawg, file.path(output_vcf_PCAWG,      paste0(chr_tag, "_PCAWG_",      suffix, ".vcf")),
            vcf_header = vcf_headers[[tag]], vcfversion = "4.2", genomereference = "hg19")
   writevcf(vcf_1kgp,  file.path(output_vcf_1kGp_PCAWG, paste0(chr_tag, "_1kGp_phase3_PCAWG_", suffix, ".vcf")),
            vcf_header = vcf_headers[[tag]], vcfversion = "4.2", genomereference = "hg19")
}

# --------------------------------------------------------------
# Main
# --------------------------------------------------------------
info       <- fread(info_path, header = TRUE, sep = "\t", data.table = FALSE)
info_males <- info[info$sexe == "male", ]

# mtDNA
vcf_MT     <- fread(path_vcf_MT, skip = "#CHROM", header = TRUE, data.table = FALSE)
vcf_MT     <- vcf_MT[grepl("VT=S(;|$)", vcf_MT$INFO), ]
vcf_MT_PCAWG <- build_vcf(vcf_MT, info$sample, "MT")
write_both_vcf(vcf_MT_PCAWG, vcf_MT, "MT", "chrMT", "snp_subset_1kgp_PCAWG_admixed_no_AMR")

# chrY (males only)
vcf_Y      <- fread(path_vcf_Y, skip = "#CHROM", header = TRUE, data.table = FALSE)
vcf_Y      <- vcf_Y[grepl("VT=SNP(;|$)", vcf_Y$INFO), ]
vcf_Y_PCAWG <- build_vcf(vcf_Y, info_males$sample, "Y")
write_both_vcf(vcf_Y_PCAWG, vcf_Y, "Y", "chrY", "snp_no_par_subset_1kgp_PCAWG_admixed_no_AMR")

# chrX
vcf_X      <- fread(path_vcf_X, skip = "#CHROM", header = TRUE, data.table = FALSE)
vcf_X      <- vcf_X[grepl("VT=SNP(;|$)", vcf_X$INFO), ]
vcf_X_PCAWG <- build_vcf(vcf_X, info$sample, "X")
write_both_vcf(vcf_X_PCAWG, vcf_X, "X", "chrX", "snp_no_par_subset_1kgp_PCAWG_admixed_no_AMR")