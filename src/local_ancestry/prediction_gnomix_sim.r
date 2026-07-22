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
SWITCH_ERROR_RATES <- c(0,
  1/1000,  # 0.001  (0.1%)
  1/500,   # 0.002  (0.2%)
  1/200,   # 0.005  (0.5%)
  1/100    # 0.01   (1%)
)


THREADS= 22
OUTPUT_DIR <- "output/local_ancestry/accuracy/prediction/1kGp_only/"
INPUT_SIM = "input/simulations"
CONFIG_PATH <- "input/gnomix"
PATH_MODEL <- "output/local_ancestry/accuracy/training/1kGp_only/" 
POP_INFO_PATH <- "rawdata/1kGp/info/1kGP.3202_samples.pop_info.txt"

# -----------------------------
# Load Data
# -----------------------------
pop_info     <- read.table(POP_INFO_PATH,    header = TRUE, sep = "\t")

ancestries_pop <- unique(na.omit(pop_info$`Super.Population`))
ancestries_pop_no_AMR <- setdiff(unique(na.omit(pop_info$`Super.Population`)), "AMR")

list_subpop <- lapply(ancestries_pop_no_AMR, function(p) {
   unique(na.omit(pop_info$Population[pop_info$`Super.Population` == p]))
})
names(list_subpop) = ancestries_pop_no_AMR
list_ancestries <- c(list(POP_no_AMR = ancestries_pop_no_AMR),list_subpop)

# -----------------------------
# Function: Prepare prediction VCF
# -----------------------------
prepare_prediction_vcf <- function(ancestries, raw_vcf_chr, output_dir_ancestry, chr) {
    cols <- fread(raw_vcf_chr, nrows = 0, skip = "#CHROM", sep = "\t", header = TRUE, data.table = FALSE)
    samples_raw <- colnames(cols)[10:length(cols)]
    check <- sapply(samples_raw, function(x) {
        x_clean <- sub("_.*$", "", x)
        parts <- strsplit(x_clean, "-")[[1]]
        all(parts %in% ancestries)
    })
                  
    samples = samples_raw[check]

    prediction_file <- if (!setequal(samples_raw, samples)) {
        file.path(output_dir_ancestry, "prediction_vcf", paste0("chr", chr, "_prediction_gnomix.vcf.gz"))
    } else {
        raw_vcf_chr
    }
  
    if (!file.exists(prediction_file) && prediction_file != raw_vcf_chr) {
    check_dir_exists_or_create(dirname(prediction_file))
    system(paste(
      "bcftools view -s",
      shQuote(paste(samples, collapse = ",")),
      shQuote(raw_vcf_chr),
      "-o", shQuote(sub("\\.gz$", "", prediction_file))
    ))
    system(paste0("gzip -f ", sub("\\.gz$", "", prediction_file)))
    }
  
  return(prediction_file)
}


# -----------------------------
# Function: Run Gnomix prediction
# -----------------------------
run_gnomix_prediction <- function(chr, seg_length, query_vcf_chr, label, output_dir_ancestry) {
    
  config_file <- file.path(CONFIG_PATH, paste0("config_", seg_length, "cM.yaml"))
    
  model_path <- file.path(PATH_MODEL, label, "training", paste0(seg_length, "cM"), paste0("chr", chr), "models",paste0("model_chm_",chr), paste0("model_chm_",chr,".pkl"))
    
  output_prediction <- file.path(output_dir_ancestry, "prediction", paste0(seg_length, "cM"), paste0("chr", chr))
  
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

lapply(names(list_ancestries), function(label){
   
    ancestries <- list_ancestries[[label]]

    lapply(SWITCH_ERROR_RATES, function(SE){
        
      
        QUERY_VCF <- file.path(INPUT_SIM, paste0("simulations_mixed_ancestry_SE_",round(SE,5),"_chrCHR.vcf.gz"))

        #mclapply(check_all_chr(QUERY_VCF), function(chr) {
        mclapply(c(22), function(chr) {

            output_dir_ancestry_SE = file.path(OUTPUT_DIR, label,paste0("SE_",round(SE,5)))
            check_dir_exists_or_create(output_dir_ancestry_SE)
            raw_vcf_chr <- gsub("chrCHR", paste0("chr", chr), QUERY_VCF)
            query_vcf_chr <- prepare_prediction_vcf(ancestries, raw_vcf_chr, output_dir_ancestry_SE, chr)

            message("Gnomix Prediction of: ", label, " - CHR",chr, " - ", SEG_LENGTH,"cM - SER = ",round(SE,5))

            run_gnomix_prediction(chr, SEG_LENGTH, query_vcf_chr, label, output_dir_ancestry_SE)
        },mc.cores=THREADS)
    })
})