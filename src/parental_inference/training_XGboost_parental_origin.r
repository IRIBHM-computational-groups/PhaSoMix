#!/usr/bin/env Rscript

# ============================================================
# Title:     XGBoost Parental Origin Model Training
# Author:    Maxime Lefebvre
# Created:   2025-11-26
# Purpose:   Train XGBoost classifiers to distinguish two ancestries
#            using germline genotypes from uniparental markers:
#              - mtDNA (maternal), chrY (paternal), chrX (maternal, males)
#            One binary model is trained per ancestry pair per marker.
#            Saved models are loaded at inference time by
#            infer_parental_origin.r via get_xgboost_prediction().
# Inputs:    - Per-marker VCF files (mtDNA, chrY, chrX)
#            - Sample info tables with 1kGp pure samples per marker
#            - Population info table (1kGP.3202_samples.pop_info.txt)
# Outputs:   output/parental_origin/training/<marker>/<pair>/xgb_model.rds
# Depends:   data.table, xgboost, utility_functions.r, XGBOOST.r
# ============================================================

suppressPackageStartupMessages(library(data.table))

source("src/utils/utility_functions.r")
source("src/utils/XGBOOST.r")

source("src/utils/XGBOOST.r")

AMR = FALSE

POP_INFO_PATH <- "rawdata/1kGp/info/1kGP.3202_samples.pop_info.txt"
OUTPUT_DIR <- "output/parental_origin/training"

INPUT_VCF_MT <- "rawdata/1kGp_PCAWG/chrMT_1kGp_phase3_PCAWG_snp_subset_1kgp_PCAWG_admixed_no_AMR.vcf.gz"
INPUT_VCF_Y <- "rawdata/1kGp_PCAWG/chrY_1kGp_phase3_PCAWG_snp_no_par_subset_1kgp_PCAWG_admixed_no_AMR.vcf.gz"
INPUT_VCF_X <- "rawdata/1kGp_PCAWG/rawdata/chrX_1kGp_phase3_PCAWG_snp_no_par_subset_1kgp_PCAWG_admixed_no_AMR.vcf.gz"

PATH_SAMPLE_INFO_MT = "rawdata/1kGp/info/samples_info_mt.txt"
PATH_SAMPLE_INFO_Y = "rawdata/1kGp/info/samples_info_y.txt"
PATH_SAMPLE_INFO_X = "rawdata/1kGp/info/samples_info_x.txt"
# -----------------------------
# Load Data
# -----------------------------

sample_info_mt_raw  = read.table(PATH_SAMPLE_INFO_MT, header=T, sep = '\t')
sample_info_y_raw  = read.table(PATH_SAMPLE_INFO_Y, header=T, sep = '\t')
sample_info_x_raw  = read.table(PATH_SAMPLE_INFO_X, header=T, sep = '\t')

if(!AMR){
    sample_info_mt_raw = sample_info_mt_raw[sample_info_mt_raw$pop != "AMR",]
    sample_info_y_raw = sample_info_y_raw[sample_info_y_raw$pop != "AMR",]
    sample_info_x_raw = sample_info_x_raw[sample_info_x_raw$pop != "AMR",]
}
  
ancestries_pop <- unique(na.omit(sample_info_mt_raw$pop))
ancestries_pop_no_AMR <- setdiff(unique(na.omit(sample_info_mt_raw$pop)), "AMR")
list_subpop <- lapply(ancestries_pop, function(p) {
   unique(na.omit(sample_info_mt_raw$subpop[sample_info_mt_raw$pop == p]))
})
names(list_subpop) <- ancestries_pop

#list_ancestries <- c(list(POP = ancestries_pop), list_subpop)
list_ancestries <- c(list(POP_no_AMR = ancestries_pop_no_AMR))


pops = unique(sample_info_mt_raw$pop)
subpops = unique(sample_info_mt_raw$subpop)

message(">>> Reading mtDNA VCF file")
vcf_mt = fread(INPUT_VCF_MT, skip = "#CHROM", header=T, sep = "\t", data.table=F)
vcf_mt <- vcf_mt[grepl("VT=S(;|$)", vcf_mt$INFO), ] # Only mono allelic SNPs

message(">>> Reading chrY VCF file")
vcf_y = fread(INPUT_VCF_Y, skip = "#CHROM", header=T, sep = "\t", data.table=F)
vcf_y <- vcf_y[grepl("VT=SNP(;|$)", vcf_y$INFO), ] # Only mono allelic SNPs 

message(">>> Reading chrX VCF file")
vcf_x = fread(INPUT_VCF_X, skip = "#CHROM", header=T, sep = "\t", data.table=F)
vcf_x <- vcf_x[grepl("VT=SNP$", vcf_x$INFO), ] # Only mono allelic SNPs
vcf_x <- vcf_x[!duplicated(vcf_x$POS), ]


# -----------------------------
# Main script
# -----------------------------


lapply(names(list_ancestries), function(cat) {
    
    message("\n")
    message("###################################################################")
    message(">>> Processing combination: ", cat)
    message("###################################################################")


    ancestries <- list_ancestries[[cat]]
    
    label <- if (all(ancestries %in% pops)) "pop" else if (all(ancestries %in% subpops)) "subpop" else NA

    # mtDNA
    message("\n  >> mtDNA")
    output_dir_mt <- file.path(OUTPUT_DIR, "mtDNA", cat)
    check_dir_exists_or_create(output_dir_mt)
    sample_file_mt  = sample_info_mt_raw[sample_info_mt_raw[[label]] %in% ancestries, c("sample", label)]     
    colnames(sample_file_mt) = c("sample", "ancestry")
    sample_file_mt$category = 0
    results_mt <- run_xgboost_pipeline(vcf_mt, sample_file_mt, output_dir_mt, "mtDNA")

    
    # chrY
    message("\n  >> chrY")
    output_dir_y <- file.path(OUTPUT_DIR, "chrY", cat)
    check_dir_exists_or_create(output_dir_y)
    sample_file_y  = sample_info_y_raw[sample_info_y_raw[[label]] %in% ancestries,c("sample", label)]     
    colnames(sample_file_y) = c("sample", "ancestry")
    sample_file_y$category = 0
    results_y <- run_xgboost_pipeline(vcf_y, sample_file_y, output_dir_y, "chrY")
    
    # chrX
    message("\n  >> chrX")
    output_dir_x <- file.path(OUTPUT_DIR, "chrX", cat)
    check_dir_exists_or_create(output_dir_x)
    sample_file_x  = sample_info_x_raw[sample_info_x_raw[[label]] %in% ancestries,c("sample", label)]     
    colnames(sample_file_x) = c("sample", "ancestry")
    sample_file_x$category = 0
    results_x <- run_xgboost_pipeline(vcf_x, sample_file_x, output_dir_x, "chrX")
}


message("\n")
message("###################################################################")
message(">>> All ancestry combinations processed!")
message("###################################################################")