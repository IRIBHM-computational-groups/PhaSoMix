#!/usr/bin/env Rscript

# ============================================================
# Title:     Phasing Switch Error Correction in Imbalanced Zones
# Author:    Maxime Lefebvre
# Created:   2025-11-26
# Purpose:   Detect and correct switch errors (SEs) in Beagle-phased
#            haplotypes for admixed cancer patients. Corrections are
#            applied only within copy-number-imbalanced zones (IZ)
#            where the allele imbalance provides a phasing signal.
#            Processing is parallelized per sample within each chromosome.
# Inputs:    - Beagle-phased VCF (per chromosome, corrected reference panel)
#            - PCAWG CNA files (ICGC and TCGA formats)
#            - Allele count files (per chromosome, per sample)
#            - VCF header file (for output VCF formatting)
# Outputs:   - Corrected per-chromosome VCF files (.vcf.gz)
#            - Summary statistics TSV (SE rates per sample/chromosome)
# Depends:   data.table, dplyr, parallel, get_BAF.r,
#            compute_SER_IS.r, writevcf.r
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(parallel)
})

source("src/utils/get_BAF.r")
source("src/utils/compute_SER_IS.r")
source("src/utils/writevcf.r")


## ============================================================
## Configuration (all paths relative to project root)
## ============================================================

input_beagle    <- "rawdata/PCAWG/rawdata/chrCHR_PCAWG_beagle5_subset_1kgp_PCAWG_snp_phased_1kGp_unrelared_IDs_high_coverage_reference_panel.vcf.gz"
cna_base_icgc   <- "/srv/home/mlef0011/rawdata/PCAWG_data/ICGC/cna/"
cna_base_tcga   <- "/srv/home/mlef0011/rawdata/PCAWG_data/TCGA/cna/"
allecount_base  <- "/mnt/iribhm/globus/secure/PCAWG_alleleCounts/"
vcf_header      <- "rawdata/vcf_headers/1kGp_hg19_header.txt"
output_vcf_dir  <- "rawdata/PCAWG/vcf_corrected/high_coverage_unrelated_1kGp/"
output_stats_path <- file.path(output_vcf_dir,"PCAWG_SE_correction_summary.txt")

## PCF and correction parameters
gamma     <- 25    # PCF penalty (higher = fewer, broader segments)
minEvents <- 5     # Minimum SNPs per PCF segment
min_run   <- 3     # Minimum consecutive LLR-negative SNPs to trigger flip

## Parallelism (set to number of available cores - 2)
N_CORES <- 20

## Create output directories
dir.create(output_vcf_dir,           showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(output_stats_path), showWarnings = FALSE, recursive = TRUE)


## ============================================================
## Per-Sample Worker Function
## ============================================================

process_sample_worker <- function(sample, chr, beagle_samples,
                                  gamma, minEvents, min_run,
                                  allecount_base, cna_base_icgc, cna_base_tcga) {

  ## Locate allele count file (per chromosome)
  allecount_file <- file.path(
    allecount_base, sample,
    gsub("chrCHR", paste0("chr", chr),
         paste0(sample, "_alleleFrequencies_chrCHR.txt.gz"))
  )

  ## Locate CNA file (try ICGC format first, then TCGA)
  cna_files <- c(
    file.path(cna_base_icgc, paste0(sample, ".consensus.20170119.somatic.cna.annotated.txt")),
    file.path(cna_base_tcga, paste0(sample, ".consensus.20170119.somatic.cna.txt"))
  )
  cna_file <- cna_files[file.exists(cna_files)][1]

  if (!file.exists(allecount_file) || length(cna_file) == 0 || is.na(cna_file)) {
    return(list(stats = NULL, pos_to_switch = c()))
  }

  ## Subset Beagle VCF to this sample
  beagle_raw <- beagle_samples[, c(colnames(beagle_samples)[1:9], sample), drop = FALSE]

  ## Compute BAF
  baf <- get_baf(sample, chr, allecount_file, beagle_raw)
  if (nrow(baf) == 0) return(list(stats = NULL, pos_to_switch = c()))

  ## Align Beagle rows to BAF positions
  beagle <- beagle_raw[match(baf$POS, beagle_raw$POS), ]

  ## Load CNA data
  cna <- tryCatch(fread(cna_file, sep = "\t", data.table = FALSE), error = function(e) NULL)
  if (is.null(cna)) return(list(stats = NULL, pos_to_switch = c()))

  ## Filter to CN-imbalanced zones on this chromosome
  iz <- subset(cna,
    chromosome == chr &
    !is.na(major_cn) & !is.na(minor_cn) &
    minor_cn != major_cn
  )
  iz <- iz[, c("chromosome", "start", "end", "major_cn", "minor_cn")]
  if (nrow(iz) == 0) return(list(stats = NULL, pos_to_switch = c()))
  if (!"total_cn" %in% colnames(iz)) iz$total_cn <- iz$major_cn + iz$minor_cn

  ## Iterate over imbalanced segments
  all_pos_to_switch <- c()
  all_pos_se        <- c()
  n_hz              <- 0

  for (idx in seq_len(nrow(iz))) {
    start     <- iz[idx, "start"]
    end       <- iz[idx, "end"]
    minor_cn  <- iz[idx, "minor_cn"]
    major_cn  <- iz[idx, "major_cn"]
    total_cn  <- iz[idx, "total_cn"]

    segment_beagle <- subset(beagle, beagle[, 1] == chr & beagle$POS >= start & beagle$POS <= end)
    segment_baf    <- subset(baf,    baf[, 1] == chr & baf$POS >= start & baf$POS <= end & baf$Good_depth > 0)

    if (nrow(segment_beagle) == 0 || nrow(segment_baf) == 0) next

    n_hz <- n_hz + sum(grepl("0\\|1|1\\|0", segment_beagle[, sample]), na.rm = TRUE)

    correction <- tryCatch(
      compute_SE_IZ(segment_baf, segment_beagle, chr,
                    minor_cn, major_cn, total_cn,
                    gamma, minEvents, min_run),
      error = function(e) NULL
    )

    if (!is.null(correction)) {
      all_pos_to_switch <- c(all_pos_to_switch, correction$pos_to_switch)
      all_pos_se        <- c(all_pos_se,        correction$pos_SE)
    }
  }

  se_rate <- if (n_hz > 0) length(all_pos_se) / n_hz * 100 else 0

  return(list(
    stats = data.frame(
      chr       = chr,
      sample    = sample,
      n_iz      = nrow(iz),
      n_switches = length(all_pos_to_switch),
      n_SE       = length(all_pos_se),
      se_rate    = round(se_rate, 2),
      n_hz_snps  = n_hz,
      stringsAsFactors = FALSE
    ),
    pos_to_switch = all_pos_to_switch
  ))
}


## ============================================================
## Main Loop (chromosomes 1–22)
## ============================================================

all_se_stats <- list()

for (chr in 1:22) {

  message(sprintf("[%s] Processing chromosome %d", Sys.time(), chr))

  beagle_file <- gsub("chrCHR", paste0("chr", chr), input_beagle)
  if (!file.exists(beagle_file)) {
    warning(sprintf("Beagle file not found for chr%d, skipping.", chr))
    next
  }

  beagle_samples <- fread(beagle_file, skip = "#CHROM", header = TRUE,
                           sep = "\t", data.table = FALSE)
  if (nrow(beagle_samples) == 0) next

  beagle_samples_corrected <- as.data.table(beagle_samples)
  samples <- colnames(beagle_samples)[10:ncol(beagle_samples)]
  message(sprintf("  Found %d samples, parallelizing over %d cores...",
                  length(samples), N_CORES))

  ## Process all samples in parallel
  se_results_chr <- mclapply(samples, function(sample) {
    process_sample_worker(sample, chr, beagle_samples,
                          gamma, minEvents, min_run,
                          allecount_base, cna_base_icgc, cna_base_tcga)
  }, mc.cores = N_CORES)
  names(se_results_chr) <- samples

  ## Apply corrections sequentially to the shared VCF data.table
  for (i in seq_along(se_results_chr)) {
    if (is.null(se_results_chr[[i]])) next
    result <- se_results_chr[[i]]
    sample <- names(se_results_chr)[i]
    if (is.null(result$stats)) next

    pos_to_switch <- result$pos_to_switch
    if (length(pos_to_switch) > 0) {
      idx_to_update <- which(beagle_samples_corrected$POS %in% pos_to_switch)
      if (length(idx_to_update) > 0) {
        beagle_samples_corrected[idx_to_update, (sample) := fifelse(
          get(sample) == "1|0", "0|1",
          fifelse(get(sample) == "0|1", "1|0", get(sample))
        )]
      }
    }
  }

  ## Collect statistics
  stats_list <- Filter(Negate(is.null), lapply(se_results_chr, `[[`, "stats"))
  if (length(stats_list) > 0) {
    se_results_chr_df <- do.call(rbind, stats_list)
    all_se_stats[[as.character(chr)]] <- se_results_chr_df
    message(sprintf("  [OK] Corrections: %d | Mean SE rate: %.2f%%",
                    sum(se_results_chr_df$n_switches),
                    mean(se_results_chr_df$se_rate, na.rm = TRUE)))
  }

  ## Write corrected chromosome VCF
  output_vcf_file <- file.path(
    output_vcf_dir,
    basename(gsub("\\.vcf\\.gz", "_corrected.vcf",
                  gsub("chrCHR", paste0("chr", chr), input_beagle)))
  )
  tryCatch({
    writevcf(as.data.frame(beagle_samples_corrected), output_vcf_file, vcf_header)
    message(sprintf("  [OK] VCF saved: %s\n", basename(output_vcf_file)))
  }, error = function(e) {
    warning(sprintf("Error writing VCF for chr%d: %s", chr, e$message))
  })

  rm(beagle_samples)
  rm(beagle_samples_corrected)

  gc()
}


## ============================================================
## Save Summary Statistics
## ============================================================

all_se_df           <- do.call(rbind, all_se_stats)
rownames(all_se_df) <- NULL

write.table(all_se_df, output_stats_path,
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

message(sprintf("\n[OK] Summary saved: %s", output_stats_path))
message(sprintf("Total samples    : %d", length(unique(all_se_df$sample))))
message(sprintf("Total corrections: %d", sum(all_se_df$n_switches)))
message(sprintf("Mean SE rate     : %.2f%%\n", mean(all_se_df$se_rate, na.rm = TRUE)))
message("Pipeline completed!")