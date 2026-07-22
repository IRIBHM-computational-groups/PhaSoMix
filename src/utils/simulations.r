#!/usr/bin/env Rscript

# ============================================================
# Title:     Mixed Ancestry Simulations with Switch Errors
# Author:    Maxime Lefebvre
# Created:   2025-03-27
# Version:   2.0
# Purpose:   Simulate mixed-ancestry individuals with
#            varying switch error rates
# Repository: https://github.com/yourusername/Phasomix
# ============================================================

## ============================================================
## Load Libraries
## ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(parallel)
})

## ============================================================
## Source Required Functions
## ============================================================
source("src/utils/writevcf.r")
source("src/utils/utility_functions.r")

## ============================================================
## Configuration
## ============================================================

# Set random seed for reproducibility
set.seed(333)

# Simulation parameters
N_MIXED_ANCESTRY <- 100           # Number of mixed-ancestry individuals per combination
SWITCH_ERROR_RATES <- c(
  1/1000,  # 0.001  (0.1%)
  1/500,   # 0.002  (0.2%)
  1/200,   # 0.005  (0.5%)
  1/100    # 0.01   (1%)
)
CHROMOSOMES <- 21:21

# File paths
OUTPUT_DIR <- "input/simulations_1kGp_only/"
#SAMPLE_INFO_PATH <- "rawdata/1kGp_PCAWG/info/sample_info_1kGp_unrelated_PCAWG_pure_split.txt"
SAMPLE_INFO_PATH <- "rawdata/1kGp/info/samples_info_1kGp_unrelated_split.txt"
POP_INFO_PATH <- "rawdata/1kGp/info/1kGP.3202_samples.pop_info.txt"
VCF_TEMPLATE <- "rawdata/1kGp_PCAWG/rawdata/chrCHR.1kGp_high_coverage_Illumina_PCAWG_beagle5_subset_1kgp_PCAWG_snp.vcf.gz"
VCF_HEADER <- "rawdata/vcf_headers/1kGp_hg19_header.txt"


## ============================================================
## Utility Functions
## ============================================================

switch_phase <- function(x) fifelse(x=="0|1", "1|0",
                             fifelse(x=="1|0", "0|1", x))

## ============================================================
## Main Simulation Functions
## ============================================================


introduce_SE <- function(vcf, samples, switch_error_rate){
    n         <- nrow(vcf)
    chr       <- unique(vcf$`#CHROM`)
    positions <- vcf$POS

    compute_switches <- function(geno) {
        idx     <- which(geno %in% c("0|1", "1|0"))
        n_het   <- length(idx)
        n_errs  <- ceiling(switch_error_rate * n_het)

        empty <- list(idx = integer(0), pos = integer(0), mask = logical(n))

        if (n_errs == 0 || n_het <= 1) return(empty)

        available <- idx[-length(idx)]
        if (length(available) == 0)      return(empty)
        n_errs <- min(n_errs, length(available))

        switch_idx <- sort(sample(available, n_errs))

        # Vectorized flip mask: each switch toggles the running parity.
        # Build a +/- step vector, cumsum it, then test odd parity.
        toggle           <- integer(n)
        toggle[switch_idx] <- 1L
        mask             <- (cumsum(toggle) %% 2L) == 1L

        list(idx = switch_idx, pos = positions[switch_idx], mask = mask)
    }

    per_sample <- lapply(samples, function(s) compute_switches(vcf[[s]]))
    names(per_sample) <- samples
                         
    vcf_SE <- vcf
    for (s in samples) {
      m <- per_sample[[s]]$mask
      if (any(m)) vcf_SE[[s]][m] <- switch_phase(vcf_SE[[s]][m])
    }
                         
    n_switches <- vapply(per_sample, function(x) length(x$idx), integer(1))

    if (sum(n_switches) == 0) {
      SE_pos <- data.frame(CHROM = character(0), INDV = character(0),
                            SER   = numeric(0),  IDX  = integer(0),
                            POS   = integer(0),  stringsAsFactors = FALSE)
    } else {
      SE_pos <- data.frame(
         CHROM = rep(chr,     sum(n_switches)),
         INDV  = rep(samples, n_switches),
         SER   = rep(switch_error_rate, sum(n_switches)),
         IDX   = unlist(lapply(per_sample, `[[`, "idx"), use.names = FALSE),
         POS   = unlist(lapply(per_sample, `[[`, "pos"), use.names = FALSE),
         stringsAsFactors = FALSE
      )
    }

    return(list(vcf_SE = vcf_SE, SE_pos = SE_pos))
                         
}


#' Create mixed ancestry individuals for one ancestry combination

create_mixed_ancestry_vcf <- function(
    vcf_test,
    ancestry_comb,
    sample_info_df,
    pop_info,
    n_samples
) {
  
  ancestries <- unlist(strsplit(ancestry_comb, '-'))
  ancestry1 <- ancestries[1]
  ancestry2 <- ancestries[2]
  
  # Determine if combination is at population or subpopulation level
  is_pop <- all(ancestries %in% unique(pop_info$Super.Population))
  is_subpop <- all(ancestries %in% unique(pop_info$Population))
  
  # Sample individuals from each ancestry
  if (is_pop) {
    samples_1 <- sample(
      sample_info_df[!is.na(sample_info_df$pop) & sample_info_df$pop == ancestry1, "sample"], 
      size = n_samples, 
      replace = FALSE
    )
    samples_2 <- sample(
      sample_info_df[!is.na(sample_info_df$pop) & sample_info_df$pop == ancestry2, "sample"], 
      size = n_samples, 
      replace = FALSE
    )
  } else if (is_subpop) {
    samples_1 <- sample(
      sample_info_df[!is.na(sample_info_df$subpop) & sample_info_df$subpop == ancestry1, "sample"], 
      size = n_samples, 
      replace = FALSE
    )
    samples_2 <- sample(
      sample_info_df[!is.na(sample_info_df$subpop) & sample_info_df$subpop == ancestry2, "sample"], 
      size = n_samples, 
      replace = FALSE
    )
  } else {
    stop(sprintf("Invalid ancestry combination: %s", ancestry_comb))
  }
  
  # Extract VCF columns for selected samples
  vcf_1 <- vcf_test[, c(colnames(vcf_test)[1:9], samples_1)]
  vcf_2 <- vcf_test[, c(colnames(vcf_test)[1:9], samples_2)]
  
  # Split phased genotypes into haplotypes
  df_1 <- cbind(
    vcf_1[, 1:9], 
    as.data.frame(t(apply(
      vcf_1[, 10:ncol(vcf_1)], 1, 
      function(x) unlist(strsplit(x, "\\|"))
    )))
  )
  
  df_2 <- cbind(
    vcf_2[, 1:9], 
    as.data.frame(t(apply(
      vcf_2[, 10:ncol(vcf_2)], 1, 
      function(x) unlist(strsplit(x, "\\|"))
    )))
  )
  
  # Shuffle haplotypes randomly
  rank_1 <- c(1:9, sample(10:ncol(df_1)))
  rank_2 <- c(1:9, sample(10:ncol(df_2)))
  
  # Create mixed ancestry individuals (one haplotype from each ancestry)
  mixed_genotypes <- lapply(10:ncol(df_1), function(i) {
    phase1 <- df_1[, rank_1[i]]
    phase2 <- df_2[, rank_2[i]]
    paste0(phase1, "|", phase2)
  })
  
  mixed_df <- do.call(cbind, mixed_genotypes)
  mixed_vcf <- cbind(df_1[, 1:9], mixed_df)
  
  # Assign sample names
  colnames(mixed_vcf)[10:ncol(mixed_vcf)] <- paste0(
    ancestry_comb, "_", 1:(n_samples * 2)
  )
  
  return(mixed_vcf)
}


#' Process one chromosome for all ancestry combinations
                                       
process_chromosome <- function(
    chr,
    combinations,
    sample_info_df,
    pop_info,
    n_ma,
    se_rates,
    output_dir
) {
  
  message("")
  message("==========================")
  message(paste0("PROCESSING CHROMOSOME", chr))
  message("==========================")
  
  ## ============================================================
  ## Load VCF data
  ## ============================================================
  
   output_simulation_vcf <- file.path(
    output_dir,
    paste0("simulations_mixed_ancestry_SE_0_chr",chr,".vcf.gz")
  )
    
  if (!file.exists(output_simulation_vcf)){
    
  message(">>> Loading VCF data...")
  
      vcf_path <- gsub("chrCHR", paste0("chr", chr), VCF_TEMPLATE)
      
      if (!file.exists(vcf_path)) {
        warning(paste0("VCF file not found: ", vcf_path))
        return(NULL)
      }

      # Filter to test samples only
      test_samples <- sample_info_df$sample
      cols_to_read <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", 
                        "FILTER", "INFO", "FORMAT", test_samples)

      vcf_test <- fread(
        vcf_path,
        header = TRUE,
        sep = "\t",
        skip = "#CHROM",
        select = cols_to_read,
        data.table = FALSE,
        showProgress = FALSE
      )

      message(paste0("  -> Loaded ", nrow(vcf_test), " variants"))
      message(paste0("  -> Samples: ",length(test_samples) ))

      ## ============================================================
      ## Simulate mixed ancestry for all combinations
      ## ============================================================

      message("")
      message(">>> Creating mixed ancestry individuals...")
      message(sprintf("  -> Ancestry combinations: %d", length(combinations)))

      sim_chr_list <- lapply(seq_along(combinations), function(i) {

        ancestry_comb <- combinations[i]
        message(sprintf("  -> Processing combination %d/%d: %s", 
                        i, length(combinations), ancestry_comb))


        ancestries <- unlist(strsplit(ancestry_comb, '-'))

        # Determine available samples
        is_pop <- all(ancestries %in% unique(pop_info$Super.Population))
        is_subpop <- all(ancestries %in% unique(pop_info$Population))

        if (is_pop) {
          n_avail_1 <- sum(!is.na(sample_info_df$pop) & sample_info_df$pop == ancestries[1])
          n_avail_2 <- sum(!is.na(sample_info_df$pop) &sample_info_df$pop == ancestries[2])
        } else if (is_subpop) {
          n_avail_1 <- sum(!is.na(sample_info_df$subpop) & sample_info_df$subpop == ancestries[1])
          n_avail_2 <- sum(!is.na(sample_info_df$subpop) & sample_info_df$subpop == ancestries[2])
        } else {
          return(NULL)
        }

        n_samples_per_group <- round(n_ma / 2)
        n_samples <- min(n_samples_per_group, n_avail_1, n_avail_2)

        if (n_samples < 1) {
          warning(sprintf("Insufficient samples for %s", ancestry_comb))
          return(NULL)
        }

        # Create mixed ancestry VCF
        mixed_vcf <- create_mixed_ancestry_vcf(
          vcf_test, ancestry_comb, sample_info_df, pop_info, n_samples
        )

        return(mixed_vcf)
      })

      # Remove NULL results
      sim_chr_list <- sim_chr_list[!sapply(sim_chr_list, is.null)]

      message(sprintf("  -> Successfully created %d combinations", length(sim_chr_list)))

      ## ============================================================
      ## Merge all combinations
      ## ============================================================

      message("")
      message(">>> Merging all ancestry combinations...")

      sim_chr_df <- Reduce(
        function(x, y) merge(x, y, by = colnames(x)[1:9], all = TRUE),
        sim_chr_list
      )

      sim_chr_df <- sim_chr_df[order(sim_chr_df$`#CHROM`, sim_chr_df$POS), ]

      message(sprintf("  -> Total simulated samples: %d", ncol(sim_chr_df) - 9))

      ## ============================================================
      ## Write baseline VCF (SE = 0)
      ## ============================================================

      message("")
      message(">>> Writing baseline VCF (no switch errors)...")

      writevcf(
        sim_chr_df,
        sub("\\.gz$", "", output_simulation_vcf),
        vcf_header = VCF_HEADER,
        genomereference = "hg19"
      )
      
  } else if (file.exists(output_simulation_vcf)){
      message(paste0("  --> ", output_simulation_vcf, " already exists. Loading file..."))
      sim_chr_df = fread(output_simulation_vcf, skip = "#CHROM", header=T, sep = "\t", data.table=F)
  }
    
  ## ============================================================
  ## Introduce switch errors at different rates
  ## ============================================================
  
  message("")
  message(">>> Introducing switch errors...")
  message(sprintf("  -> Processing %d error rates ...", length(se_rates))) 
    
    samples <- colnames(sim_chr_df)[10:ncol(sim_chr_df)]

    lapply(se_rates, function(se_rate) {

      se_vcf_path <- file.path(output_dir, paste0("simulations_mixed_ancestry_SE_",round(se_rate,5),"_chr",chr,".vcf.gz"))
      se_pos_path <- file.path(output_dir, paste0("SE_positions_simulations_mixed_ancestry_",round(se_rate,5),"_chr",chr,".tsv"))

      if (file.exists(se_vcf_path) && file.exists(se_pos_path)) return(NULL)
     
      message(paste0(" -> Switch error rate: ", round(se_rate,5)))
        
      SE_info = introduce_SE(sim_chr_df,samples, se_rate)
        
      write.table(SE_info$SE_pos, se_pos_path, col.names=T, row.names=F, quote=F, sep = "\t")
      writevcf(SE_info$vcf_SE, sub("\\.gz$", "", se_vcf_path), vcf_header = VCF_HEADER, genomereference = "hg19")
        
    })
  
  message("")
  message(sprintf(">>> Chromosome %d completed!", chr))
  
  return(NULL)
}


## ============================================================
## Main Execution
## ============================================================

run_mixed_ancestry_simulations <- function() {
  
  message("====================================================")
  message("MIXED ANCESTRY SIMULATIONS WITH SWITCH ERRORS")
  message("====================================================")
  message(sprintf("Started at: %s", Sys.time()))
  message("")
  
  ## ============================================================
  ## Load data and generate combinations
  ## ============================================================
  
  message(">>> Loading configuration...")
  message(sprintf("  -> Number of mixed-ancestry individuals: %d", N_MIXED_ANCESTRY))
  message(sprintf("  -> Switch error rates: %d", length(SWITCH_ERROR_RATES)))
  message(sprintf("  -> Chromosomes to process: %d", length(CHROMOSOMES)))
  
  # Create output directory
  dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
  
  # Load population information
  message("")
  message(">>> Loading population data...")
  pop_info <- read.table(POP_INFO_PATH, sep = '\t', header = TRUE)
  message(sprintf("  -> Populations: %d", length(unique(pop_info$Super.Population))))
  message(sprintf("  -> Subpopulations: %d", length(unique(pop_info$Population))))
  
  # Load sample information (test samples only)
  samples_info = read.table(SAMPLE_INFO_PATH, header=T, sep = "\t")

  sample_info_df <- as.data.frame(subset(samples_info, category == 1))
  message(sprintf("  -> Test samples: %d", nrow(sample_info_df)))
  
  # Generate ancestry combinations
  combinations <- generate_ancestry_combinations(pop_info)
  
  ## ============================================================
  ## Process each chromosome
  ## ============================================================
  
  message("")
  message("====================================================")
  message("STARTING CHROMOSOME PROCESSING")
  message("====================================================")
  
  results <- lapply(CHROMOSOMES, function(chr) {
    process_chromosome(
      chr = chr,
      combinations = combinations,
      sample_info_df = sample_info_df,
      pop_info = pop_info,
      n_ma = N_MIXED_ANCESTRY,
      se_rates = SWITCH_ERROR_RATES,
      output_dir = OUTPUT_DIR
    )
  })
  
  ## ============================================================
  ## Summary
  ## ============================================================
  
  message("")
  message("====================================================")
  message("SIMULATIONS COMPLETE")
  message("====================================================")
  message(sprintf("Completed at: %s", Sys.time()))
  message("")
  message("Summary:")
  message(sprintf("  Chromosomes processed: %d", length(CHROMOSOMES)))
  message(sprintf("  Ancestry combinations: %d", length(combinations)))
  message(sprintf("  Switch error rates: %d", length(SWITCH_ERROR_RATES)))
  message(sprintf("  Output directory: %s", OUTPUT_DIR))
  message("")
  message("Output files per chromosome:")
  message("  - simulations_mixed_ancestry_SE_0_chrN.vcf.gz")
  for (se in SWITCH_ERROR_RATES) {
    message(sprintf("  - simulations_mixed_ancestry_SE_%.5f_chrN.vcf.gz", se))
    message(sprintf("  - SE_pos_simulations_mixed_ancestry_SE_%.5f_chrN.tsv", se))
  }
  message(paste(rep("=", 70), collapse = ""))
  print(warnings())
  invisible(NULL)
}

## ============================================================
## Execute
## ============================================================


run_mixed_ancestry_simulations()