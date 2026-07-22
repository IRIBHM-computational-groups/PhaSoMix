#!/usr/bin/env Rscript

# ============================================================
# Title:     B-Allele Frequency (BAF) Computation
# Author:    Maxime Lefebvre
# Created:   2025-11-26
# Purpose:   Compute per-SNP BAF for each phased haplotype using
#            allele counts and Beagle phasing. BAF is used by
#            compute_SER_IS to detect switch errors in CN-imbalanced
#            genomic zones.
# Inputs:    - Beagle phased VCF (data.frame, already loaded)
#            - Allele count file (tab-separated, per chromosome)
# Outputs:   BAF data frame (optionally saved and/or plotted)
# Depends:   data.table, ggplot2, utility_functions.r
# ============================================================

suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
})

source("src/utils/utility_functions.r")


## ============================================================
## Phasing Helper
## ============================================================

#' Extract phase1/phase2 allele labels from a Beagle VCF data frame
#'
#' Keeps only heterozygous sites (0|1 or 1|0), then maps the
#' integer genotype (0=REF, 1=ALT) to the actual nucleotide.
#'
#' @param x Data frame: Beagle VCF with at least columns
#'           #CHR, POS, REF, ALT, and one sample column (col 10).
#' @return Data frame with columns: chr, pos, phase1, phase2.
readPhases <- function(x) {
  # Keep heterozygous phased SNPs only
  x[, 10] <- sub(":.*", "", x[, 10])
  x       <- x[x[, 10] %in% c("0|1", "1|0"), ]

  phase1 <- gsub("(.*)\\|(.*)", "\\1", x[, 10])
  phase2 <- gsub("(.*)\\|(.*)", "\\2", x[, 10])

  # Map 0 -> REF allele, 1 -> ALT allele
  phases1 <- x[, "REF"]
  phases1[phase1 == "1"] <- x[phase1 == "1", "ALT"]
  phases2 <- x[, "REF"]
  phases2[phase2 == "1"] <- x[phase2 == "1", "ALT"]

  data.frame(
    chr    = x[, 1],
    pos    = x[, 2],
    phase1 = phases1,
    phase2 = phases2
  )
}


## ============================================================
## Main BAF Function
## ============================================================

#' Compute B-Allele Frequency for each phased haplotype
#'
#' For each heterozygous SNP, BAF_phase1 = count(allele in phase1) / total_depth.
#' The function supports on-disk caching (rewrite = FALSE) and optional
#' QC plots saved alongside the BAF file.
#'
#' @param sample          Character. Sample identifier.
#' @param chr             Integer or character. Chromosome number (e.g. 1).
#' @param allele_count_path Character. Path to the allele count file for this
#'                        chromosome (tab-separated with columns #CHR, POS,
#'                        Count_A/C/G/T, Good_depth).
#' @param beagle          Data frame. Beagle phased VCF already loaded in memory
#'                        (rows = SNPs, sample column at position 10).
#' @param output_path     Character. Directory where BAF files are saved.
#' @param save            Logical. Write BAF table to disk (default FALSE).
#' @param rewrite         Logical. Overwrite existing BAF file (default TRUE).
#' @param plot            Logical. Save a BAF scatter plot (default FALSE).
#' @return Data frame with columns: #CHR, POS, phase1, phase2,
#'         Count_A/C/G/T, Good_depth, BAF_phase1, BAF_phase2.
get_baf <- function(sample, chr, allele_count_path, beagle,
                    output_path = ".", save = TRUE,
                    rewrite = FALSE, plot = TRUE) {

  output_baf_dir <- file.path(output_path, "baf")
  check_dir_exists_or_create(output_baf_dir)

  output_baf <- file.path(output_baf_dir,
                          paste0("baf_chr", as.character(chr), "_", sample, ".txt"))

  ## Return cached result if available
  if (file.exists(output_baf) && !rewrite) {
    baf <- fread(output_baf, sep = "\t", header = TRUE, data.table = FALSE)
    return(baf)
  }

  message("  >> BAF computing [CHR", chr, "]...")

  ## Load allele counts
  allele_count <- fread(allele_count_path, header = TRUE, sep = "\t", data.table = FALSE)

  ## Extract phased alleles from Beagle VCF
  phases <- readPhases(beagle)
  colnames(phases)[1:2] <- colnames(allele_count)[1:2]
  phases <- phases[order(phases$POS, decreasing = FALSE), ]

  ## Merge on chromosome + position
  baf <- merge(phases, allele_count, by = c("#CHR", "POS"))

  ## Compute BAF for phase 1
  baf$BAF_phase1 <- 0
  for (nuc in c("A", "C", "G", "T")) {
    baf$BAF_phase1[baf$phase1 == nuc] <- baf[[paste0("Count_", nuc)]][baf$phase1 == nuc]
  }
  baf$BAF_phase1 <- baf$BAF_phase1 / baf$Good_depth

  ## Compute BAF for phase 2
  baf$BAF_phase2 <- 0
  for (nuc in c("A", "C", "G", "T")) {
    baf$BAF_phase2[baf$phase2 == nuc] <- baf[[paste0("Count_", nuc)]][baf$phase2 == nuc]
  }
  baf$BAF_phase2 <- baf$BAF_phase2 / baf$Good_depth

  baf[is.na(baf)] <- 0
  baf <- baf[order(baf$POS, decreasing = FALSE), ]

  ## Optionally save to disk
  if (save) {
    write.table(baf, output_baf, sep = "\t",
                col.names = TRUE, row.names = FALSE, quote = FALSE)
  }

  ## Optionally generate QC plot
  if (plot) {
    baf_plot <- ggplot(data = baf, aes(x = POS, y = BAF_phase1)) +
      geom_point(size = 0.05, fill = "grey") +
      geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
      scale_x_continuous(labels = function(x) x / 1e6, name = "Position (Mb)") +
      labs(title = paste0(sample, " - chr", chr, " BAF (phase 1)")) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5))

    ggsave(filename = gsub("\\.txt$", ".png", output_baf),
           plot = baf_plot, width = 6, height = 4, dpi = 300)
  }

  return(baf)
}