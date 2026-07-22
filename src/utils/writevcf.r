#!/usr/bin/env Rscript

# ============================================================
# Title:     VCF File Writer
# Author:    Maxime Lefebvre
# Created:   2025-11-26
# Purpose:   Write a data frame as a properly formatted VCF file
#            with header, sorted variants, and gzip compression.
#            Used after phase-correction to produce output VCFs.
# Inputs:    - Data frame with standard VCF columns
#            - Optional path to an existing header file
# Outputs:   Gzip-compressed VCF file (.vcf.gz)
# Depends:   data.table, utility_functions.r
# ============================================================

source("src/utils/utility_functions.r")

suppressWarnings(suppressMessages(library(data.table)))


## ============================================================
## Main VCF Writing Function
## ============================================================

writevcf <- function(input,
                     vcf_path,
                     vcf_header    = FALSE,
                     vcfversion    = "4.2",
                     genomereference = "hg19") {

  message(">>> Writing VCF file...")

  ## ---- Validate required columns ----
  required_columns <- c("#CHROM", "POS", "ID", "REF", "ALT",
                        "QUAL", "FILTER", "INFO", "FORMAT")
  input <- check_required_columns(input, required_columns)

  ## ---- Sort by chromosome and position ----
  if (inherits(input, "data.table")) {
    if (!is.numeric(input$POS)) input[, POS := as.numeric(POS)]
    setorder(input, `#CHROM`, POS)
  } else if (inherits(input, "data.frame")) {
    if (!is.numeric(input$POS)) input$POS <- as.numeric(input$POS)
    input <- input[order(input$`#CHROM`, input$POS), ]
  } else {
    stop("Input must be a data.frame or data.table", call. = FALSE)
  }

  ## ---- Write header ----
  if (isFALSE(vcf_header)) {
    message("  -> Generating default VCF header")
    header_lines <- c(
      sprintf("##fileformat=VCFv%s",   vcfversion),
      sprintf("##fileDate=%s",          format(Sys.Date(), "%Y%m%d")),
      sprintf("##reference=%s",         genomereference),
      '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
    )
    cat(paste(header_lines, collapse = "\n"), "\n", file = vcf_path)

  } else if (is.character(vcf_header) && file.exists(vcf_header)) {
    message(sprintf("  -> Copying header from: %s", vcf_header))
    system_result <- system(paste("cat", shQuote(vcf_header), ">", shQuote(vcf_path)))
    if (system_result != 0) stop("Failed to copy VCF header", call. = FALSE)

  } else {
    stop("vcf_header must be FALSE or a path to an existing header file",
         call. = FALSE)
  }

  ## ---- Write variant data ----
  message("  -> Writing variant data...")
  suppressWarnings(
    write.table(input, file = vcf_path,
                sep = "\t", col.names = TRUE, row.names = FALSE,
                quote = FALSE, append = TRUE)
  )

  ## ---- Compress ----
  message("  -> Compressing with gzip...")
  gzip_result <- system(sprintf("gzip -f %s", shQuote(vcf_path)))
  output_path <- if (gzip_result == 0) {
    paste0(vcf_path, ".gz")
  } else {
    warning("gzip compression failed, keeping uncompressed file")
    vcf_path
  }

  ## ---- Verify output ----
  if (!file.exists(output_path)) {
    stop("VCF file was not created successfully", call. = FALSE)
  }
  message(sprintf("  -> VCF successfully written: %s", output_path))
  message(">>> VCF writing complete!")

  invisible(output_path)
}