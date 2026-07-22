#!/usr/bin/env Rscript

# ============================================================
# Title:     Regional SNV Quantification
# Author:    Maxime Lefebvre
# Created:   2025-12-01
# Purpose:   Compute CN-corrected SNV counts in specific genomic
#            regions (exons, introns, intergenic, genes, CpG,
#            imprinted loci) for one sample. Each region is
#            identified by overlapping SNV positions with BED
#            or GFF annotation files using GenomicRanges.
# Inputs:    - SNV list (VCF-style TSV)
#            - SNV quantification file (output of SNV_correction.r)
#            - Reference genome (BSgenome / DNAStringSet)
#            - CNA file (passed through to get_stats)
# Outputs:   Per-region TSV files + get_stats() results
# Depends:   GenomicRanges, Biostrings, dplyr, data.table, get_stats.r
# ============================================================

suppressPackageStartupMessages({
    library(GenomicRanges)
    library(dplyr)
    library(Biostrings)
    library(data.table)
})

source("src/utils/get_stats.r")


## ============================================================
## Annotation File Paths (relative to project root)
## ============================================================
## These can be overridden by passing explicit paths to
## local_quantification().

input_exon_path        <- "rawdata/reference/exon_regions_hg19_sorted.bed"
input_intron_path      <- "rawdata/reference/intron_regions_hg19.bed"
input_intergenic_path  <- "rawdata/reference/intergenic_regions_hg19_sorted.bed"
input_genes_path       <- "rawdata/reference/gencode.v19.annotation_sorted.gff"
input_imprinted_path   <- "rawdata/reference/Imprintome_ICRs_hg19.txt"


## ============================================================
## CpG Context Filter
## ============================================================

get_C2TatCpG <- function(sample, input_SNV, input_SNV_quant, genome, output_path) {

  SNV      <- fread(input_SNV,       sep = "\t", header = TRUE, data.table = FALSE)
  SNV = SNV[SNV$FILTER == ".",]
  SNV_quant <- fread(input_SNV_quant, sep = "\t", header = TRUE, data.table = FALSE)

  ## Retrieve flanking dinucleotides from reference
  SNV_quant$binucleotide_minus <- mapply(function(chr, pos) {
    as.character(unlist(subseq(genome[[paste0("chr", chr)]], start = pos - 1, end = pos)))
  }, SNV_quant$chr_snv, SNV_quant$pos_snv)

  SNV_quant$binucleotide_plus <- mapply(function(chr, pos) {
    as.character(unlist(subseq(genome[[paste0("chr", chr)]], start = pos, end = pos + 1)))
  }, SNV_quant$chr_snv, SNV_quant$pos_snv)

  data <- merge(SNV[, c("#CHROM", "POS", "REF", "ALT")], SNV_quant,
                by.x = c("#CHROM", "POS"), by.y = c("chr_snv", "pos_snv"))
  data$transition <- paste0(data$REF, ">", data$ALT)

  ## C>T at CpG (forward strand) or G>A at CpG (reverse strand)
  C2TatCpG <- data[
    (data$binucleotide_plus  == "CG" & data$transition == "C>T") |
    (data$binucleotide_minus == "CG" & data$transition == "G>A"),
  ]

  output_path_C2TatCpG <- file.path(output_path,
                                     paste0("SNV_quantification_C2TatCpG_", sample, ".txt"))
  write.table(C2TatCpG, output_path_C2TatCpG,
              quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

  return(output_path_C2TatCpG)
}


## ============================================================
## Generic BED Overlap Filter
## ============================================================

get_region_bed <- function(bed_file, input_SNV_quant, output_path) {

  colnames(bed_file) <- c("chr", "start", "end")

  SNVs <- fread(input_SNV_quant, sep = "\t", data.table = FALSE)
  SNVs$chr_snv <- paste0("chr", as.character(SNVs$chr_snv))

  snvs_gr <- GRanges(
    seqnames = SNVs$chr_snv,
    ranges   = IRanges(start = SNVs$pos_snv, end = SNVs$pos_snv)
  )
  bedfile_gr <- GRanges(
    seqnames = bed_file$chr,
    ranges   = IRanges(start = bed_file$start, end = bed_file$end)
  )

  overlaps      <- findOverlaps(snvs_gr, bedfile_gr)
  filtered_snvs <- SNVs[queryHits(overlaps), ] %>% distinct()
  filtered_snvs$chr_snv <- gsub("chr", "", filtered_snvs$chr_snv)

  write.table(filtered_snvs, output_path,
              quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
}


## ============================================================
## Main Regional Quantification Function
## ============================================================

local_quantification <- function(sample,
                                 input_SNV,
                                 input_SNV_quant,
                                 genome,
                                 input_CNA,          
                                 output_path,
                                 input_exon       = input_exon_path,
                                 input_intron     = input_intron_path,
                                 input_intergenic = input_intergenic_path,
                                 input_genes      = input_genes_path,
                                 input_imprinted  = input_imprinted_path) {

  ## ---- CpG context ----
  C2TatCpG_path <- get_C2TatCpG(sample, input_SNV, input_SNV_quant, genome, output_path)

  ## ---- Exons ----
  bed_file <- read.table(input_exon)
  output_path_exon <- file.path(output_path, paste0("SNV_quantification_exon_", sample, ".txt"))
  get_region_bed(bed_file, input_SNV_quant, output_path_exon)

  ## ---- Introns ----
  bed_file <- read.table(input_intron)
  output_path_intron <- file.path(output_path, paste0("SNV_quantification_intron_", sample, ".txt"))
  get_region_bed(bed_file, input_SNV_quant, output_path_intron)

  ## ---- Intergenic ----
  bed_file <- read.table(input_intergenic)
  output_path_intergenic <- file.path(output_path, paste0("SNV_quantification_intergenic_", sample, ".txt"))
  get_region_bed(bed_file, input_SNV_quant, output_path_intergenic)

  ## ---- Genes (from GFF) ----
  bed_file <- read.table(input_genes)
  bed_file <- bed_file[bed_file$V3 == "gene", c("V1", "V4", "V5")]
  output_path_genes <- file.path(output_path, paste0("SNV_quantification_genes_", sample, ".txt"))
  get_region_bed(bed_file, input_SNV_quant, output_path_genes)

  ## ---- Imprinted loci ----
  bed_file <- read.table(input_imprinted, header = TRUE)
  bed_file <- bed_file[
    !is.na(bed_file$chr) &
    !is.na(bed_file$"Parental.Origin.of.methylation") &
    bed_file$"Parental.Origin.of.methylation" != "S",
  ]
  output_path_imprinted <- file.path(output_path, paste0("SNV_quantification_imprinted_", sample, ".txt"))
  get_region_bed(bed_file, input_SNV_quant, output_path_imprinted)

  ## ---- Compute statistics for each region ----
  stats_CpG       <- get_stats(sample, C2TatCpG_path, input_CNA, output_path, "CpG")
  stats_exon      <- get_stats(sample, output_path_exon, input_CNA, output_path, "exon")
  stats_intron    <- get_stats(sample, output_path_intron, input_CNA, output_path, "intron")
  stats_intergenic <- get_stats(sample, output_path_intergenic, input_CNA, output_path, "intergenic")
  stats_genes     <- get_stats(sample, output_path_genes, input_CNA, output_path, "genes")
  stats_imprinted <- get_stats(sample, output_path_imprinted, input_CNA, output_path, "imprinted")

  return(list(
    stats_CpG        = stats_CpG,
    stats_exon       = stats_exon,
    stats_intron     = stats_intron,
    stats_intergenic = stats_intergenic,
    stats_genes      = stats_genes,
    stats_imprinted  = stats_imprinted
  ))
}