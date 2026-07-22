#!/usr/bin/env Rscript

# ============================================================
# Title:     Compute ancestry informativeness scores per SNP
# Author:    Maxime Lefebvre
# Created:   19/05/2025
# Purpose:   Parse 1kG high-coverage VCFs hg19 to calculate
#            the ancestry informativeness of each SNP based on 
#            allele frequencies (AF) across populations.
# Requires:  R packages: data.table, stringr, parallel
# Output:    One score file per chromosome + global file listing all outputs
# ============================================================

# ---------------------- Load libraries ----------------------
library(data.table)
library(stringr)
library(parallel)
library(optparse)

# ---------------------- Parse command line arguments ----------------------
option_list <- list(
  make_option(c("-i", "--input"), type = "character", help = "Input VCF path with 'chrX' placeholder"),
  make_option(c("-o", "--output"), type = "character", help = "Output path with 'chrX' placeholder"),
  make_option(c("-t", "--threads"), type = "integer", default = 10, help = "Number of threads [default %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

input_template <- opt$input
output_template <- opt$output
threads <- opt$threads

# ---------------------- Define utility functions ----------------------

compute_score <- function(col) {
  best <- head(sort(col, decreasing = TRUE), 2)
  best_idx <- head(order(col, decreasing = TRUE), 2)
  ancestries <- c("EAS_AF", "AFR_AF", "AMR_AF", "EUR_AF", "SAS_AF")
  best_ancestry <- ancestries[best_idx]

  max1 <- exp(best[1]) / sum(exp(col))
  max2 <- exp(best[2]) / sum(exp(col))
  res <- max1 - max2

  return(list(Result = res, best_ancestry = best_ancestry[1]))
}

get_score <- function(input_vcf, output){
  vcf_raw <- fread(input_vcf, skip = "#CHROM", sep = "\t", header = TRUE, select = 1:9)
  vcf_snp <- vcf_raw[grepl("VT=SNP$", vcf_raw$INFO), ]

  AFs <- c("EAS_AF", "AFR_AF", "AMR_AF", "EUR_AF", "SAS_AF")
  AF_values <- lapply(AFs, function(AF)
    as.numeric(str_extract(vcf_snp$INFO, paste0("(?<=", AF, "=)[0-9.]+")))
  )
  AF_values_df <- as.data.frame(do.call(rbind, AF_values))

  AF_score <- apply(AF_values_df, 2, compute_score)
  AF_score_df <- as.data.frame(do.call(rbind, lapply(AF_score, unlist)))

  AF_score_df <- cbind(vcf_snp, AF_score_df)
  colnames(AF_score_df)[10:11] <- c("score", "top_ancestry")
  AF_score_df$score <- as.numeric(AF_score_df$score)

  AF_score_df <- AF_score_df[order(AF_score_df$score, decreasing = TRUE), ]
  write.table(AF_score_df, output, col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

  return(AF_score_df)
}

# ---------------------- Process all chromosomes ----------------------

all_chr_results <- lapply(1:22, function(chr) {

  cat("\n>>> Processing CHROMOSOME", chr, "...\n")

  input_vcf <- sub("chrX", paste0("chr", chr), input_template)
  output <- sub("chrX", paste0("chr", chr), output_template)

  if (!file.exists(input_vcf)) {
    message("File not found for chr", chr, ": ", input_vcf, " — skipping.")
    return(NULL)
  }

  score_vcf <- tryCatch({
    get_score(input_vcf, output)
  }, error = function(e) {
    message("Error processing chr", chr, ": ", conditionMessage(e))
    return(NULL)
  })

  return(score_vcf)

})

# ---------------------- Save combined results ----------------------

all_chr_df <- do.call(rbind, all_chr_results[!sapply(all_chr_results, is.null)])
summary_output <- sub("chrX", "all_chr", output_template)
write.table(all_chr_df, summary_output, col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
