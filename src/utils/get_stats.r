#!/usr/bin/env Rscript

# ============================================================
# Title:     SNV Count Statistics per Ancestry and CN Zone
# Author:    Maxime Lefebvre
# Created:   2025-12-01
# Purpose:   Compute CN-corrected SNV counts for each ancestry
#            copy (copy1/copy2) and parental origin (mother/father),
#            stratified by copy-number zone:
#              - no_IZ    : balanced regions
#              - IZ_no_LOH: imbalanced, minor_cn > 0
#              - IZ_LOH   : loss of heterozygosity (minor_cn == 0)
# Inputs:    - Quantification file (output of SNV_correction.r)
# Outputs:   Named vector written to disk as TSV
# Depends:   data.table
# ============================================================

suppressPackageStartupMessages({
    library(data.table)
})

get_stats <- function(sample, quantification_SNV_path, input_CNA,
                      output_path, label) {

  data <- fread(quantification_SNV_path, sep = "\t", header = TRUE)

  ## ---- Helper: sum corrected counts for a subset ----
  sumcorr <- function(filt) round(sum(as.numeric(data$correction[filt]), na.rm = TRUE))

  ## ---- Total SNVs ----
  nbr_mut_hap1     <- sumcorr(data$haplotype_SNV == "hapA")
  nbr_mut_hap2     <- sumcorr(data$haplotype_SNV == "hapB")

  ## ---- Total SNVs ----
  nbr_mut_hap1_no_loh     <- sumcorr(data$haplotype_SNV == "hapA" & data$type != "IZ_LOH")
  nbr_mut_hap2_no_loh    <- sumcorr(data$haplotype_SNV == "hapB" & data$type != "IZ_LOH")
    
  ## ---- Non-imbalanced zones (balanced CN) ----
  nbr_mut_hap1_no_iz     <- sumcorr(data$type == "no_IZ" & data$haplotype_SNV == "hapA")
  nbr_mut_hap2_no_iz     <- sumcorr(data$type == "no_IZ" & data$haplotype_SNV == "hapB")

  ## ---- Imbalanced zones (any IZ) ----
  iz_mask <- data$type %in% c("IZ_LOH", "IZ_no_LOH")
  nbr_mut_hap1_iz     <- sumcorr(iz_mask & data$haplotype_SNV == "hapA")
  nbr_mut_hap2_iz     <- sumcorr(iz_mask & data$haplotype_SNV == "hapB")

  ## ---- IZ with LOH ----
  nbr_mut_hap1_iz_LOH     <- sumcorr(data$type == "IZ_LOH" & data$haplotype_SNV == "hapA")
  nbr_mut_hap2_iz_LOH     <- sumcorr(data$type == "IZ_LOH" & data$haplotype_SNV == "hapB")

  ## ---- IZ without LOH ----
  nbr_mut_hap1_iz_no_loh     <- sumcorr(data$type == "IZ_no_LOH" & data$haplotype_SNV == "hapA")
  nbr_mut_hap2_iz_no_loh     <- sumcorr(data$type == "IZ_no_LOH" & data$haplotype_SNV == "hapB")

  ## ---- Assemble result vector ----
  res <- c(
    sample            = sample,

    SNVs_hap1        = nbr_mut_hap1,
    SNVs_hap2        = nbr_mut_hap2,
      
    SNVs_hap1_no_LOH     = nbr_mut_hap1_no_loh,
    SNVs_hap2_no_LOH     = nbr_mut_hap2_no_loh,

    SNVs_no_IZ_hap1  = nbr_mut_hap1_no_iz,
    SNVs_no_IZ_hap2  = nbr_mut_hap2_no_iz,

    SNVs_IZ_no_LOH_hap1  = nbr_mut_hap1_iz_no_loh,
    SNVs_IZ_no_LOH_hap2  = nbr_mut_hap2_iz_no_loh
  )

  write.table(res,
              file.path(output_path, paste0("stats_SNV_quantification_", label, ".txt")),
              col.names = FALSE, row.names = TRUE, sep = "\t", quote = FALSE)

  return(res)
}