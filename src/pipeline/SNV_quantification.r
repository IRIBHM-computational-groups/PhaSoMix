#!/usr/bin/env Rscript

# ============================================================
# SNV Quantification Pipeline
# Purpose: Phase SNVs, infer ancestry, quantify mutations
# ============================================================

suppressPackageStartupMessages({
    library(optparse)
    library(Biostrings)
    library(parallel)
})

# ---- Command Line Arguments ----
option_list <- list(
  make_option(c("-i", "--path_info_ma"), type = "character",
    help = "TSV with columns: sample, ancestry, sexe, input_bam, input_allelecount, input_SNV, input_CNA, input_SNV_ccf"),
  
  make_option(c("-o", "--output_path"), type = "character", default = getwd(),
    help = "Output directory [default: current directory]"),
  
  make_option(c("-v", "--vcf_ma"), type = "character",
    help = "VCF path with 'CHR' placeholder (e.g., vcf_chrCHR.vcf.gz)"),
    
  make_option(c("-s", "--vcf_ma_sex"), type = "character",
    help = "VCF path with 'CHR' placeholder (e.g., vcf_chrCHR.vcf.gz)"),
  
  make_option(c("-r", "--ref_genome_path"), type = "character",
    help = "Reference genome FASTA path"),
  
  make_option(c("-g", "--gnomix_pred_path"), type = "character",
    help = "Gnomix prediction directory"),
  
  make_option(c("-x", "--xgboost_model_path"), type = "character",
    default = "output/parental_origin/training/",
    help = "XGBoost model directory [default: output/parental_origin/training/]"),
  
  make_option(c("-a", "--accuracy"), type = "character",
    default = "output/accuracy_PhaSomix.txt",
    help = "Accuracy file path [default: output/accuracy_PhaSomix.txt]"),
      
  make_option(c("-t", "--threads"), type = "integer",
    default = 10,
    help = "Number of threads  [default:10]")
)

parser <- OptionParser(option_list = option_list,
  description = "SNV Quantification Pipeline")
args <- parse_args(parser)

# ---- Load Required Functions ----
source("src/utils/utility_functions.r")
source("src/pipeline/SNV_phasing.r")
source("src/pipeline/SNV_ancestry.r")
source("src/pipeline/SNV_correction.r")
source("src/utils/get_stats.r")
source("src/utils/local_quantification.r")
source("src/utils/get_circos_plot.r")

# ---- Load Input Data ----
info_ma <- read.table(args$path_info_ma, header = TRUE, sep = "\t")

# Validate required columns
info_ma <- check_required_columns(info_ma, 
  c('sample', 'project', 'hapA', 'hapB', 'ancestry', 'sexe', 'input_bam', 'input_allelecount', 
    'input_SNV', 'input_CNA', 'input_SNV_ccf'))

# ---- Set Parameters ----
output_path <- args$output_path 
ref_genome_path <- args$ref_genome_path
input_pred_gnomix <- args$gnomix_pred_path 
input_model_xgboost <- args$xgboost_model_path 
input_vcf <- args$vcf_ma 
input_vcf_sex <- args$vcf_ma_sex 
THREADS <- args$threads

# Load reference genome (chr1-22 only)
genome <- readDNAStringSet(ref_genome_path)[1:22]
names(genome) <- paste0("chr", 1:22)

samples <- unique(info_ma$sample)

# ---- Process Each Sample ----
results <- lapply(samples, function(sample_ID) {
  
  tryCatch({
    
    message("\n========================================")
    message("Processing sample: ", sample_ID)
    message("========================================\n")
    
    # Extract sample-specific info
    label <- unlist(strsplit(info_ma[info_ma$sample == sample_ID, 'ancestry'], "-"))
    hapA = sort(unlist(strsplit(gsub("[{}]", "", info_ma[info_ma$sample == sample_ID, 'hapA']), ",")))
    hapB = sort(unlist(strsplit(gsub("[{}]", "", info_ma[info_ma$sample == sample_ID, 'hapB']), ",")))
    sexe <- info_ma[info_ma$sample == sample_ID, "sexe"]
    input_bam <- info_ma[info_ma$sample == sample_ID, "input_bam"]
    input_allelecount <- info_ma[info_ma$sample == sample_ID, "input_allelecount"]
    input_SNV <- info_ma[info_ma$sample == sample_ID, "input_SNV"]
    input_CNA <- info_ma[info_ma$sample == sample_ID, "input_CNA"]
    input_SNV_ccf <- info_ma[info_ma$sample == sample_ID, "input_SNV_ccf"]
    
    # Check all files exist
    files_to_check <- c(input_bam, input_SNV, input_CNA, input_SNV_ccf)
    missing_files <- files_to_check[!file.exists(files_to_check)]
    
    if (length(missing_files) > 0) {
      message("SKIPPED: Missing files - ", paste(missing_files, collapse = ", "))
      return(NULL)
    }
    
    # Create output directory
    output_path_sample <- file.path(output_path, sample_ID)
    check_dir_exists_or_create(output_path_sample)
    
    # Step 1: Phase SNVs using nearby heterozygous SNPs
    message(">>> STEP 1: SNV PHASING")
    phased_SNV_path <- phase_SNV(1:22, sample_ID, input_bam, input_vcf, 
                                  output_path_sample, input_SNV, ref_genome_path,rewrite=FALSE)
    
    # Step 2: Infer ancestry for each phased SNV
    message("\n>>> STEP 2: ANCESTRY INFERENCE")

    ancestry_SNV_path <- get_ancestry(sample_ID, hapA, hapB, 1:22,
                                      phased_SNV_path, input_pred_gnomix, 
                                      input_vcf, output_path_sample, verbose = FALSE, rewrite=FALSE)
      
    
    # Step 3: Quantify SNVs with copy number correction
    message("\n>>> STEP 4: QUANTIFICATION & CN CORRECTION")
    quantification_SNV_path <- get_quantification(sample_ID, ancestry_SNV_path, 
                                                  input_SNV_ccf, input_vcf, 
                                                  input_allelecount, output_path_sample, rewrite=FALSE)
    
    # Step 4: Generate circos plot
    message("\n>>> STEP 5: CIRCOS PLOT")
    get_circos_plot(sample_ID, ancestry_SNV_path, input_SNV, input_CNA, 
                   input_pred_gnomix, output_path_sample, CNV_MAX = 5, ref_genome = "hg19", rewrite=FALSE)
          
      
     
    # Step 5: Compute statistics for different genomic regions
    message("\n>>> STEP 6: REGIONAL STATISTICS")
    stats <- c(
      list(stats_global = get_stats(sample_ID, quantification_SNV_path, 
                                    input_CNA, output_path_sample, "global")),
      local_quantification(sample_ID, input_SNV, quantification_SNV_path, 
                          genome, input_CNA, output_path_sample)
    
    )
    message("\n“ Sample ", sample_ID, " completed successfully\n")
    return(stats)
    }, error = function(e) {
    message("\n— ERROR for sample ", sample_ID, ": ", e$message, "\n")
    return(NULL)
  })
})

message("\n Pipeline completed successfully!")