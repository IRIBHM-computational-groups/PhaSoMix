#!/usr/bin/env Rscript

# ============================================================
# Title:     Training Gnomix (multi-class)
# Author:    Maxime Lefebvre
# Created:   2025/11/26
# ============================================================

source("src/utils/utility_functions.r")
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(parallel)))


# -----------------------------
# Global Parameters
# -----------------------------
set.seed(333)
SEG_LENGTHS <- c(0.1)

OUTPUT_DIR       <- "output/local_ancestry/accuracy/training/1kGp_only/"
#SAMPLE_INFO_PATH <- "rawdata/1kGp_PCAWG/info/sample_info_1kGp_unrelated_PCAWG_pure_split.txt"
SAMPLE_INFO_PATH <- "rawdata/1kGp/info/samples_info_1kGp_unrelated_split.txt"
CONFIG_PATH      <- "input/gnomix/"
QUERY_VCF        <- "None"
RAW_VCF          <- "rawdata/1kGp_PCAWG/rawdata/chrCHR.1kGp_high_coverage_Illumina_PCAWG_beagle5_subset_1kgp_PCAWG_snp.vcf.gz"
GENETIC_MAP      <- "rawdata/1kGp/mapfiles/plink.ALLchr.GRCh37.map"

# If TRUE, downsample each ancestry to the size of the smallest one
# so the training set is class-balanced.
BALANCE_CLASSES  <- TRUE

THREADS <- 2


# -----------------------------
# Load Data
# -----------------------------
samples_info <- read.table(SAMPLE_INFO_PATH, header = TRUE, sep = "\t")

# NB: the column names below must match samples_info; adjust if needed.
ancestries_pop <- unique(na.omit(samples_info$pop))
ancestries_pop_no_AMR <- setdiff(unique(na.omit(samples_info$pop)), "AMR")
list_subpop <- lapply(ancestries_pop_no_AMR, function(p) {
   unique(na.omit(samples_info$subpop[samples_info$pop == p]))
})
names(list_subpop) <- ancestries_pop_no_AMR


list_ancestries <- c(list(POP_no_AMR = ancestries_pop_no_AMR), list_subpop)
list_ancestries = list_subpop
# Cache chromosome list so we do not re-read the VCF header per loop.
ALL_CHRS <- check_all_chr(RAW_VCF)
ALL_CHRS=22

# -----------------------------
# Function: Create sample map
# -----------------------------
create_sample_map <- function(ancestries, samples_info, balance = TRUE) {

   is_pop    <- all(ancestries %in% unique(na.omit(samples_info$pop)))
   is_subpop <- all(ancestries %in% unique(na.omit(samples_info$subpop)))

   if (is_pop) {
      sample_map <- subset(samples_info,
                           pop %in% ancestries & category == 0)
      sample_map <- sample_map[, c("sample", "pop")]
   } else if (is_subpop) {
      sample_map <- subset(samples_info,
                           subpop %in% ancestries & category == 0)
      sample_map <- sample_map[, c("sample", "subpop")]
   } else {
      return(NULL)
   }

   colnames(sample_map) <- c("#Sample", "Panel")

   sample_map <- sample_map[!is.na(sample_map$Panel) &
                            sample_map$Panel != "", ]

   if (isTRUE(balance)) {
      counts <- table(sample_map$Panel)
      n_min  <- min(counts)
      message(sprintf("  -> Class sizes before balancing: %s",
                      paste(sprintf("%s=%d", names(counts), counts),
                            collapse = ", ")))
      message(sprintf("  -> Balancing to %d samples per class", n_min))

      sample_map <- do.call(rbind, lapply(ancestries, function(a) {
         sub <- sample_map[sample_map$Panel == a, ]
         sub[sample(seq_len(nrow(sub)), n_min), ]
      }))
   }

   missing_anc <- setdiff(ancestries, unique(sample_map$Panel))
   if (length(missing_anc) > 0) {
      warning("Missing ancestries in sample_map: ",
              paste(missing_anc, collapse = ", "))
      return(NULL)
   }

   sample_map
}


# -----------------------------
# Function: Prepare training VCF
# -----------------------------
prepare_training_vcf <- function(sample_map, raw_vcf_chr, output_dir, chr) {
   cols <- fread(raw_vcf_chr, nrows = 0, skip = "#CHROM", sep = "\t",
                 header = TRUE, data.table = FALSE)
   samples_raw <- as.vector(unlist(cols[10:length(cols)]))

   training_file <- if (!setequal(sample_map$`#Sample`, samples_raw)) {
      file.path(output_dir, "training_vcf",
                paste0("chr", chr, "_training_gnomix.vcf.gz"))
   } else {
      raw_vcf_chr
   }

   if (!file.exists(training_file) && training_file != raw_vcf_chr) {
      check_dir_exists_or_create(dirname(training_file))
      system(paste(
         "bcftools view -s",
         shQuote(paste(sample_map$`#Sample`, collapse = ",")),
         shQuote(raw_vcf_chr),
         "-o", shQuote(sub("\\.gz$", "", training_file))
      ))
      system(paste0("gzip -f ", sub("\\.gz$", "", training_file)))
   }

   training_file
}


# -----------------------------
# Function: Run Gnomix training
# -----------------------------
run_gnomix_training <- function(chr, seg_length, training_file,
                                output_sample_file, output_dir) {

   config_file <- file.path(CONFIG_PATH,
                            paste0("config_", seg_length, "cM.yaml"))

   output_training <- file.path(output_dir, "training",
                                paste0(seg_length, "cM"),
                                paste0("chr", chr))
   check_dir_exists_or_create(output_training)

   cmd <- paste(
      "python3",
      "software/gnomix/gnomix.py",
      shQuote(QUERY_VCF),
      shQuote(output_training),
      shQuote(chr),
      "False",
      shQuote(GENETIC_MAP),
      shQuote(training_file),
      shQuote(output_sample_file),
      shQuote(config_file),
      "2>&1 | tee",
      shQuote(file.path(output_training,
                        paste0("chr", chr, "_", seg_length,
                               "cM_training_log.txt")))
   )

   cat("  -> Running command:\n\t", cmd, "\n")
   system(cmd)
   cat("  -> Chromosome", chr, "completed.\n")
}


# -----------------------------
# Main: iterate over training sets (POP, AFR, EAS, EUR, SAS, ...)
# -----------------------------
lapply(names(list_ancestries), function(label) {

   ancestries <- list_ancestries[[label]]

   message("\n=== Training run: ", label,
           " (", length(ancestries), " classes: ",
           paste(ancestries, collapse = ", "), ") ===")

   sample_map <- create_sample_map(ancestries, samples_info,
                                   balance = BALANCE_CLASSES)
   if (is.null(sample_map)) {
      warning("Could not build sample_map for ", label, "; skipping.")
      return(NULL)
   }
   message(sprintf("  -> sample_map: %d samples across %d classes",
                   nrow(sample_map), length(ancestries)))

   output_dir <- file.path(OUTPUT_DIR, label)
   check_dir_exists_or_create(output_dir)

   output_sample_file <- file.path(
      output_dir,
      paste0("sample_map_gnomix_training_", label, ".smap")
   )
   write.table(sample_map, output_sample_file,
               col.names = TRUE, row.names = FALSE,
               sep = "\t", quote = FALSE)

   mclapply(ALL_CHRS, function(chr) {
      raw_vcf_chr   <- gsub("chrCHR", paste0("chr", chr), RAW_VCF)
      training_file <- prepare_training_vcf(sample_map, raw_vcf_chr,
                                            output_dir, chr)

      lapply(SEG_LENGTHS, function(seg_length) {
         run_gnomix_training(chr, seg_length, training_file,
                             output_sample_file, output_dir)
      })
   }, mc.cores = THREADS)
})

message("\nAll training runs done.")