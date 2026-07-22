#!/usr/bin/env Rscript

# ============================================================
# Title:     Discriminant SNP Extraction
# Author:    Maxime Lefebvre
# Created:   2025-11-26
# Purpose:   Select the most informative (discriminant) SNPs for
#            ADMIXTURE-based global ancestry inference. SNPs are
#            filtered by MAF and ancestry-specific scores, then
#            extracted from per-chromosome VCF files.
# Inputs:    - Per-chromosome VCF files (raw_vcf template with 'CHR')
#            - SNP score file (precomputed ancestry informativeness)
#            - Sample information table
# Outputs:   Merged VCF (.vcf.gz) with selected discriminant SNPs
# Depends:   data.table, dplyr, stringr, writevcf.r, utility_functions.r
# ============================================================

source("src/utils/writevcf.r")
source("src/utils/utility_functions.r")

suppressWarnings(suppressMessages({
  library(stringr)
  library(dplyr)
  library(data.table)
}))


## ============================================================
## Allele Frequency Helper
## ============================================================

#' Extract allele frequency from a VCF INFO field string
#' @param info_string Character vector of VCF INFO fields.
#' @return Numeric vector of AF values.
extract_af <- function(info_string) {
  as.numeric(str_extract(info_string, "(?<=AF=)[0-9.]+"))
}


## ============================================================
## SNP Filtering
## ============================================================

#' Filter discriminant SNPs by MAF and ancestry-specific score
#'
#' Selects the top-scoring SNPs for each ancestry group with equal
#' representation, after applying a MAF threshold.
#'
#' @param score_file        Path to the precomputed SNP score file.
#' @param ancestries        Character vector of ancestry labels.
#' @param maf               Numeric. Minimum minor allele frequency threshold.
#' @param n_snps            Integer. Total target number of SNPs to select.
#' @param output_directory  Character. Directory for writing output files.
#' @return Data frame of selected discriminant SNPs.
filter_informative_snps <- function(score_file, ancestries, maf,
                                    n_snps, output_directory = ".") {

  output_file <- file.path(
    output_directory,
    sprintf("discriminant_SNPs_score_best_%d_MAF_tsh_%.3f.txt", n_snps, maf)
  )

  ## Return cached result
  if (file.exists(output_file)) {
    discr_snps <- fread(output_file, header = TRUE, sep = "\t")
    message(sprintf(">>> Loaded cached discriminant SNPs (%d SNPs)", nrow(discr_snps)))
    return(discr_snps)
  }

  message(">>> Filtering informative SNPs...")
  check_dir_exists_or_create(output_directory)

  score_df <- fread(score_file, header = TRUE, sep = "\t", data.table = FALSE)

  ## Compute MAF from INFO field
  message("  -> Extracting allele frequencies...")
  score_df$AF  <- extract_af(score_df$INFO)
  score_df$MAF <- pmin(score_df$AF, 1 - score_df$AF)

  n_total <- nrow(score_df)
  message(sprintf("  -> Total SNPs: %d | Filtering MAF > %.3f...", n_total, maf))

  score_df_filtered <- subset(score_df, MAF > maf)
  n_after_maf <- nrow(score_df_filtered)
  message(sprintf("  -> SNPs after MAF filter: %d (%.1f%%)",
                  n_after_maf, 100 * n_after_maf / n_total))

  ## Check that all requested ancestries are present in score file
  ancestries_AF     <- paste0(ancestries, "_AF")
  available         <- unique(score_df_filtered$top_ancestry)
  missing_ancestries <- setdiff(ancestries_AF, available)
  if (length(missing_ancestries) > 0) {
    stop(sprintf("Ancestries not found in score file: %s",
                 paste(sub("_AF$", "", missing_ancestries), collapse = ", ")),
         call. = FALSE)
  }

  score_df_filtered <- subset(score_df_filtered, top_ancestry %in% ancestries_AF)

  ## Select top SNPs per ancestry (equal representation)
  n_ancestries      <- length(unique(score_df_filtered$top_ancestry))
  snps_per_ancestry <- floor(n_snps / n_ancestries)
  message(sprintf("  -> Selecting %d SNPs per ancestry...", snps_per_ancestry))

  discr_snps <- score_df_filtered %>%
    group_by(top_ancestry) %>%
    arrange(desc(score)) %>%
    slice_head(n = snps_per_ancestry) %>%
    ungroup() %>%
    arrange(desc(score))

  message("  -> SNPs per ancestry:")
  for (anc in names(table(discr_snps$top_ancestry))) {
    message(sprintf("     - %s: %d",
                    sub("_AF$", "", anc), table(discr_snps$top_ancestry)[anc]))
  }

  write.table(discr_snps, output_file,
              col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
  message(sprintf("  -> Saved: %s", output_file))
  message(sprintf(">>> Filtering complete! Selected %d SNPs", nrow(discr_snps)))

  return(discr_snps)
}


## ============================================================
## VCF Extraction
## ============================================================

#' Extract discriminant SNPs from per-chromosome VCF files
#'
#' Full pipeline:
#'   1. Filter informative SNPs (MAF + ancestry score)
#'   2. Extract those SNPs from each VCF chromosome file
#'   3. Merge across chromosomes and write a single VCF
#'
#' @param vcf_files        Character vector of per-chromosome VCF paths.
#' @param score_file       Path to SNP score file.
#' @param samples_info     Data frame with column 'sample'.
#' @param ancestries       Character vector of ancestry labels.
#' @param maf              Numeric MAF threshold.
#' @param n_snps           Integer target number of SNPs.
#' @param output_directory Character output directory.
#' @return Path to the output VCF.gz file.
extract_discriminant_SNPs <- function(vcf_files, score_file, samples_info,
                                      ancestries, maf, n_snps,
                                      output_directory) {

  message(strrep("=", 50))
  message("DISCRIMINANT SNP EXTRACTION PIPELINE")
  message(strrep("=", 50))
  message(sprintf("Started at: %s", Sys.time()))

  check_dir_exists_or_create(output_directory)

  ## STEP 1: Filter discriminant SNPs
  message("\n[STEP 1/2] Filtering informative SNPs...")
  discriminant_snps <- filter_informative_snps(
    score_file, ancestries, maf, n_snps, output_directory
  )

  ## STEP 2: Extract from VCF files
  output_vcf_path <- file.path(
    output_directory,
    sprintf("discriminant_SNPs_%d_%.3f.vcf.gz", n_snps, maf)
  )

  if (file.exists(output_vcf_path)) {
    message(sprintf("[STEP 2/2] Cached output found: %s", output_vcf_path))
    return(output_vcf_path)
  }

  message("\n[STEP 2/2] Extracting SNP genotypes from VCFs...")

  discr_snp_list <- lapply(seq_along(vcf_files), function(i) {
    vcf_path <- vcf_files[i]
    message(sprintf("  -> File %d/%d: %s", i, length(vcf_files), basename(vcf_path)))

    vcf <- fread(vcf_path, skip = "#CHROM", sep = "\t",
                 header = TRUE, data.table = FALSE, showProgress = FALSE)

    vcf_samples      <- colnames(vcf)[10:ncol(vcf)]
    samples_to_keep  <- intersect(vcf_samples, samples_info$sample)
    missing_samples  <- setdiff(samples_info$sample, vcf_samples)

    if (length(missing_samples) > 0) {
      warning(sprintf("  -> Missing %d samples from %s",
                      length(missing_samples), basename(vcf_path)))
    }
    message(sprintf("     Keeping %d samples", length(samples_to_keep)))

    vcf <- vcf[, c(colnames(vcf)[1:9], samples_to_keep)]

    ## Filter to discriminant SNP positions
    snps_chr <- merge(
      vcf,
      discriminant_snps[, c("#CHROM", "POS", "REF", "ALT")],
      by = c("#CHROM", "POS", "REF", "ALT")
    )
    snps_chr <- snps_chr[order(snps_chr$POS), ]
    message(sprintf("     Extracted %d SNPs", nrow(snps_chr)))
    return(snps_chr)
  })

  ## Combine across chromosomes
  common_cols  <- Reduce(intersect, lapply(discr_snp_list, colnames))
  discr_snp_df <- do.call(rbind, lapply(discr_snp_list,
                                        function(df) df[, common_cols, drop = FALSE]))

  ## Remove columns with NA values
  cols_with_na <- names(discr_snp_df)[colSums(is.na(discr_snp_df)) > 0]
  if (length(cols_with_na) > 0) {
    message(sprintf("  -> Removing %d columns with NA", length(cols_with_na)))
    discr_snp_df <- discr_snp_df[, !(names(discr_snp_df) %in% cols_with_na)]
  }

  ## Write merged VCF
  message("\n>>> Writing merged VCF...")
  writevcf(discr_snp_df, sub("\\.gz$", "", output_vcf_path))

  ## Summary
  message("\n", strrep("=", 50))
  message("EXTRACTION COMPLETE")
  message(sprintf("  VCFs processed : %d", length(vcf_files)))
  message(sprintf("  SNPs selected  : %d", nrow(discriminant_snps)))
  message(sprintf("  Samples output : %d", ncol(discr_snp_df) - 9))
  message(sprintf("  Output VCF     : %s", output_vcf_path))
  message(strrep("=", 50))

  return(output_vcf_path)
}