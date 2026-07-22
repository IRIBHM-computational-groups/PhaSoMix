#!/usr/bin/env Rscript

# ============================================================
# Title:     Prediction Gnomix Paired
# Author:    Maxime Lefebvre
# Created:   2025/11/26
# Purpose:   Predict local ancestry from pre-trained Gnomix models 
#            for each pairwaise of test samples
# ============================================================

# Load utility functions
source("src/utils/utility_functions.r")
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(parallel)))

# -----------------------------
# Global Parameters
# -----------------------------
SEG_LENGTH <- 0.1

THREADS= 22
OUTPUT_DIR <- "output/local_ancestry/prediction/high_coverage/IS_no_corrected/"
CONFIG_PATH <- "input/gnomix/"
SAMPLE_INFO_PATH <- "rawdata/1kGp_PCAWG/info/sample_info_1kGp_unrelated_PCAWG_pure_PCAWG_admixed.txt"
PATH_MODEL <- "output/local_ancestry/training/" 
QUERY_VCF <- "/srv/home/mlef0011/Phasomix/rawdata/1kGp_PCAWG/rawdata/chrCHR.1kGp_high_coverage_Illumina_PCAWG_beagle5_subset_1kgp_PCAWG_snp.vcf.gz"

# -----------------------------
# Load Data
# -----------------------------
samples_info <- read.table(SAMPLE_INFO_PATH, header = TRUE, sep = "\t")

# NB: the column names below must match samples_info; adjust if needed.
ancestries_pop <- unique(na.omit(samples_info$pop))
ancestries_pop_no_AMR <- setdiff(unique(na.omit(samples_info$pop)), "AMR")
list_subpop <- lapply(ancestries_pop, function(p) {
   unique(na.omit(samples_info$subpop[samples_info$pop == p]))
})
names(list_subpop) <- ancestries_pop

#list_ancestries <- c(list(POP = ancestries_pop), list_subpop)
#list_ancestries <- c(list(POP_no_AMR_no_balanced = ancestries_pop))
list_ancestries <- c(list(POP_no_AMR = ancestries_pop_no_AMR))

# Cache chromosome list so we do not re-read the VCF header per loop.
ALL_CHRS <- check_all_chr(QUERY_VCF)


# -----------------------------
# Function: Prepare prediction VCF
# -----------------------------

prepare_prediction_vcf <- function(samples_info, raw_vcf_chr, output_dir_ancestry, chr) {
  cols <- fread(raw_vcf_chr, nrows = 0, skip = "#CHROM", sep = "\t", header = TRUE, data.table = FALSE)
  samples_raw <- colnames(cols)[10:length(cols)]
  
  samples = samples_info$sample[samples_info$category == 1]
    
  prediction_file <- if (!setequal(samples, samples_raw)) {
    file.path(output_dir_ancestry, "prediction_vcf", paste0("chr", chr, "_prediction_gnomix.vcf.gz"))
  } else {
    raw_vcf_chr
  }
  
  if (!file.exists(prediction_file) && prediction_file != raw_vcf_chr) {
    check_dir_exists_or_create(dirname(prediction_file))
    system(paste(
      "bcftools view --force-samples -s",
      shQuote(paste(samples, collapse = ",")),
      shQuote(raw_vcf_chr),
      "-o", shQuote(sub("\\.gz$", "", prediction_file))
    ))
    system(paste0("bgzip -f ", sub("\\.gz$", "", prediction_file)))
  }
  
  return(prediction_file)
}

# -----------------------------
# Function: Run Gnomix prediction
# -----------------------------
run_gnomix_prediction <- function(chr, seg_length, query_vcf_chr, label, output_dir) {
    
  config_file <- file.path(CONFIG_PATH, paste0("config_", seg_length, "cM.yaml"))
    
  model_path <- file.path(PATH_MODEL, label, "training", paste0(seg_length, "cM"), paste0("chr", chr), "models",paste0("model_chm_",chr), paste0("model_chm_",chr,".pkl"))
    
  output_prediction <- file.path(output_dir, "prediction", paste0(seg_length, "cM"), paste0("chr", chr))
  
  check_dir_exists_or_create(output_prediction)
  
  cmd <- paste(
    "python3",
    "software/gnomix/gnomix.py",
    shQuote(query_vcf_chr),
    shQuote(output_prediction),
    shQuote(chr),
    "False",
    shQuote(model_path),
    shQuote(config_file),
    "2>&1 | tee",
    shQuote(file.path(output_prediction, paste0("chr", chr, "_", seg_length, "cM_prediction_log.txt")))
  )
  
  message("  -> Running command:\n\t", cmd, "\n")
  system(cmd)
  message("  -> Chromosome ", chr, " completed.\n")
}

# -----------------------------
# Main Processing
# -----------------------------

lapply(names(list_ancestries), function(label) {
   
   ancestries <- list_ancestries[[label]]
        
    mclapply(ALL_CHRS, function(chr) {

        output_dir = file.path(OUTPUT_DIR, label)
        check_dir_exists_or_create(output_dir)
        raw_vcf_chr <- gsub("chrCHR", paste0("chr", chr), QUERY_VCF)
            
        query_vcf_chr <- prepare_prediction_vcf(samples_info, raw_vcf_chr, output_dir, chr)

        message("Gnomix Prediction of: ", label, " - CHR",chr, " - ", SEG_LENGTH,"cM")

        run_gnomix_prediction(chr, SEG_LENGTH, query_vcf_chr, label, output_dir)
    },mc.cores=THREADS)
})
