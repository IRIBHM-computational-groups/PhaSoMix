#!/usr/bin/env Rscript

# ============================================================
# Title:     ADMIXTURE Global Ancestry Inference Wrapper
# Author:    Maxime Lefebvre
# Created:   2025-11-26
# Purpose:   Run supervised global ancestry estimation with ADMIXTURE
#            (Alexander et al. 2009) on a cohort of samples, at both
#            population level (e.g. AFR, EUR) and sub-population level
#            (e.g. YRI, CEU). The pipeline has three steps:
#              1. Convert multi-sample VCF to PLINK binary format
#                 (.bed/.bim/.fam) via PLINK, retaining only samples
#                 listed in sample_info.
#              2. Generate per-analysis .pop label files: training
#                 samples (category = 0) receive their known ancestry
#                 label; test samples (category = 1) receive "-" so
#                 that ADMIXTURE infers their proportions.
#              3. Run ADMIXTURE in supervised mode. Parse the Q-matrix
#                 (ancestry proportion per sample) and write a
#                 formatted TSV with sample IDs and ancestry columns.
#            The main entry point is get_admixture_prediction(), which
#            optionally first runs extract_discriminant_SNPs() to reduce
#            the marker set to the most informative sites.
# Inputs:    - Multi-sample VCF (discriminant SNPs, or raw with CHR placeholder)
#            - sample_info TSV: sample, pop, subpop, category (0=train, 1=test)
#            - score_file: per-SNP ancestry discrimination scores
# Outputs:   - <output_dir>/pop/<prefix>.<K>.Q (ADMIXTURE Q matrix)
#            - <output_dir>/pop/<prefix>_admixture_prediction_pop.txt
#            - <output_dir>/subpop/<prefix>_admixture_prediction_subpop.txt
# Depends:   utility_functions.r, PLINK (system), ADMIXTURE (system),
#            extract_discriminant_SNPs.r
# ============================================================

source("src/utils/utility_functions.r")
source("src/utils/extract_discriminant_SNPs.r")


## ============================================================
## Default Executable Paths (override via function arguments)
## ============================================================

DEFAULT_PLINK_PATH     <- "software/plink"
DEFAULT_ADMIXTURE_PATH <- "software/admixture"


## ============================================================
## Step 1: VCF → PLINK Binary Conversion
## ============================================================

#' Convert a multi-sample VCF to PLINK binary format
#'
#' Calls PLINK --vcf --make-bed to produce .bed/.bim/.fam files.
#' Only samples listed in sample_info are retained (--keep).
#' Two copies of the PLINK files are created: one under output_path/pop/
#' for population-level ADMIXTURE and one under output_path/subpop/
#' for sub-population-level ADMIXTURE.
#'
#' @param input_vcf   Character. Path to the input VCF (may be gzipped).
#' @param output_path Character. Base output directory.
#' @param sample_info Data frame. Must have columns: sample, and at least
#'                    one of pop or subpop.
#' @param exec_plink  Character. Path to PLINK executable.
#' @return Named list with elements pop and/or subpop, each a character
#'         string giving the PLINK file prefix (without extension).
make_bed <- function(input_vcf, output_path, sample_info,
                     exec_plink = DEFAULT_PLINK_PATH) {

    message(">>> Converting VCF to PLINK binary format...")

    ## -- Validation --
    if (!file.exists(input_vcf))  stop("Input VCF not found: ", input_vcf,  call. = FALSE)
    if (!file.exists(exec_plink)) stop("PLINK executable not found: ", exec_plink, call. = FALSE)

    has_pop    <- "pop"    %in% colnames(sample_info)
    has_subpop <- "subpop" %in% colnames(sample_info)
    if (!has_pop && !has_subpop)
        stop("sample_info must have a 'pop' and/or 'subpop' column", call. = FALSE)

    ## -- Build temporary conversion directory --
    base_name        <- sub("\\.vcf.*$", "", basename(input_vcf))
    temp_dir         <- file.path(output_path, "temp_bed")
    check_dir_exists_or_create(temp_dir)
    temp_output_file <- file.path(temp_dir, base_name)

    message(sprintf("  -> Input VCF: %s", input_vcf))
    message(sprintf("  -> Samples to retain: %d", nrow(sample_info)))

    ## -- Write sample keep file (FID IID format required by PLINK) --
    keep_file <- file.path(temp_dir, "keep_samples.txt")
    write.table(data.frame(FID = sample_info$sample, IID = sample_info$sample),
                keep_file, quote = FALSE, row.names = FALSE, col.names = FALSE)

    ## -- Run PLINK --
    plink_cmd <- paste(exec_plink,
                       "--vcf",    shQuote(input_vcf),
                       "--keep",   shQuote(keep_file),
                       "--make-bed",
                       "--double-id ",
                       "--out",    shQuote(temp_output_file),
                       "--silent")
    message("  -> Running: ", plink_cmd)
    if (system(plink_cmd) != 0) stop("PLINK conversion failed", call. = FALSE)

    plink_ext <- c(".bed", ".bim", ".fam")
    if (!all(file.exists(paste0(temp_output_file, plink_ext))))
        stop("PLINK did not produce expected output files", call. = FALSE)
    message("  -> PLINK files created.")

    ## -- Copy to pop/ and subpop/ analysis directories --
    paths <- list()
    for (level in intersect(c("pop", "subpop"), colnames(sample_info))) {
        level_dir  <- file.path(output_path, level)
        check_dir_exists_or_create(level_dir)
        paths[[level]] <- file.path(level_dir, base_name)
        file.copy(paste0(temp_output_file, plink_ext),
                  paste0(paths[[level]], plink_ext), overwrite = TRUE)
        message(sprintf("  -> Copied to %s/", level))
    }

    system(paste("rm -rf", shQuote(temp_dir)))
    message(">>> PLINK conversion complete.")
    return(paths)
}


## ============================================================
## Step 2: Generate .pop Label Files for Supervised ADMIXTURE
## ============================================================

#' Create ADMIXTURE .pop files for supervised mode
#'
#' Reads the PLINK .fam file to get the sample order, then writes a
#' .pop file where training samples (category == 0) receive their
#' known ancestry label and test samples (category == 1) receive "-"
#' (tells ADMIXTURE to infer proportions for those samples).
#'
#' @param input_plink Named list (from make_bed). Keys: pop, subpop.
#'                    Values: PLINK file prefixes (without extension).
#' @param sample_info Data frame with columns: sample, pop, subpop, category.
#' @return Named list with elements output_pop and output_subpop
#'         (paths to written .pop files, or NA if level not used).
make_pop_file <- function(input_plink, sample_info) {

    message(">>> Generating ADMIXTURE population files...")

    output_pop    <- NA
    output_subpop <- NA

    ## Internal helper: write one .pop file for a given analysis level
    write_pop <- function(plink_prefix, label_col, level_name) {
        fam      <- read.table(paste0(plink_prefix, ".fam"), header = FALSE,
                                stringsAsFactors = FALSE)
        si_ord   <- sample_info[match(fam$V1, sample_info$sample), ]
        labels   <- ifelse(si_ord$category == 0, si_ord[[label_col]], "-")
        out_path <- paste0(plink_prefix, ".pop")
        writeLines(as.character(labels), out_path, useBytes = TRUE)

        n_train  <- sum(si_ord$category == 0)
        n_test   <- sum(si_ord$category == 1)
        n_labels <- length(unique(si_ord[[label_col]][si_ord$category == 0]))
        message(sprintf("  -> %s: %d training | %d test | %d %s classes → %s",
                         level_name, n_train, n_test, n_labels, level_name, out_path))
        return(out_path)
    }

    if (!is.null(input_plink$pop))    output_pop    <- write_pop(input_plink$pop,    "pop",    "pop")
    if (!is.null(input_plink$subpop)) output_subpop <- write_pop(input_plink$subpop, "subpop", "subpop")

    message(">>> Population files ready.")
    return(list(output_pop = output_pop, output_subpop = output_subpop))
}


## ============================================================
## Step 3: Run Supervised ADMIXTURE
## ============================================================

#' Execute ADMIXTURE in supervised mode and parse Q-matrix output
#'
#' Determines K from the number of unique non-"-" labels in the .pop
#' file, then runs ADMIXTURE --supervised. ADMIXTURE writes its output
#' to the current working directory, so the Q and P files are renamed
#' to their intended locations.
#' The Q-matrix is merged with sample IDs and ancestry labels, then
#' written as a tab-separated TSV.
#'
#' @param input_plink Named list. PLINK file prefixes (pop and/or subpop).
#' @param pop_file    Named list (from make_pop_file). Keys: output_pop,
#'                    output_subpop. Values: .pop file paths.
#' @param output_path Character. Base output directory.
#' @param admx_exec   Character. Path to ADMIXTURE executable.
#' @return Named list with elements output_pop and output_subpop
#'         (paths to written prediction TSVs).
run_ADMIXTURE <- function(input_plink, pop_file, output_path,
                           admx_exec = DEFAULT_ADMIXTURE_PATH) {

    message(">>> Running supervised ADMIXTURE...")
    if (!file.exists(admx_exec)) stop("ADMIXTURE executable not found: ", admx_exec, call. = FALSE)

    ## Internal helper: run ADMIXTURE for one analysis level
    run_one_level <- function(pop_vec, sample_ids, plink_prefix, level_name) {

        output_dir <- file.path(output_path, level_name)
        check_dir_exists_or_create(output_dir)

        ## K = number of reference populations (excluding "-")
        pops    <- unique(pop_vec[pop_vec != "-"])
        K       <- length(pops)
        message(sprintf("  -> %s: K = %d (%s)", level_name, K, paste(pops, collapse = ", ")))

        ## ADMIXTURE writes output to the current working directory
        admx_cmd <- paste(admx_exec, "--supervised", shQuote(paste0(plink_prefix, ".bed")), K)
        message("  -> Command: ", admx_cmd)
        if (system(admx_cmd, intern = FALSE) != 0)
            stop(sprintf("ADMIXTURE failed for %s", level_name), call. = FALSE)

        ## Rename ADMIXTURE output files from cwd to their final locations
        tmp_q <- file.path(getwd(), paste0(basename(plink_prefix), ".", K, ".Q"))
        tmp_p <- file.path(getwd(), paste0(basename(plink_prefix), ".", K, ".P"))
        fin_q <- paste0(plink_prefix, ".", K, ".Q")
        fin_p <- paste0(plink_prefix, ".", K, ".P")
        file.rename(tmp_q, fin_q)
        file.rename(tmp_p, fin_p)
        if (!file.exists(fin_q)) stop("Expected Q file not found: ", fin_q, call. = FALSE)

        ## Parse Q-matrix: each row = one sample, each column = one ancestry proportion
        q_mat   <- read.table(fin_q, header = FALSE)
        results <- cbind(data.frame(sample_ID = sample_ids, ancestry = pop_vec), q_mat)
        colnames(results) <- c("sample_ID", "ancestry", as.character(pops))

        out_file <- file.path(output_dir,
                               sprintf("%s_admixture_prediction_%s.txt",
                                       basename(plink_prefix), level_name))
        write.table(results, out_file, col.names = TRUE, row.names = FALSE,
                    sep = "\t", quote = FALSE)
        message(sprintf("  -> Written: %s  (%d samples)", out_file, nrow(results)))
        return(out_file)
    }

    output_pop    <- NA
    output_subpop <- NA

    if (!is.na(pop_file$output_pop)) {
        fam        <- read.table(paste0(input_plink$pop, ".fam"), header = FALSE)
        pop_vec    <- scan(pop_file$output_pop, what = character(), quiet = TRUE)
        output_pop <- run_one_level(pop_vec, fam$V1, input_plink$pop, "pop")
    }
    if (!is.na(pop_file$output_subpop)) {
        fam           <- read.table(paste0(input_plink$subpop, ".fam"), header = FALSE)
        subpop_vec    <- scan(pop_file$output_subpop, what = character(), quiet = TRUE)
        output_subpop <- run_one_level(subpop_vec, fam$V1, input_plink$subpop, "subpop")
    }

    message(">>> ADMIXTURE complete.")
    return(list(output_pop = output_pop, output_subpop = output_subpop))
}


## ============================================================
## Full ADMIXTURE Pipeline (Steps 1–3)
## ============================================================

#' Run the complete ADMIXTURE pipeline from VCF to predictions
#'
#' Convenience wrapper that calls make_bed() → make_pop_file() →
#' run_ADMIXTURE() in sequence and prints a summary.
#'
#' @param input_vcf   Character. Path to input VCF.
#' @param sample_info Data frame. Sample metadata (sample, pop, subpop, category).
#' @param output_dir  Character. Root output directory.
#' @param exec_plink  Character. Path to PLINK executable.
#' @param admx_exec   Character. Path to ADMIXTURE executable.
#' @return Named list with output_pop and output_subpop prediction file paths.
ADMIXTURE <- function(input_vcf, sample_info, output_dir,
                       exec_plink = DEFAULT_PLINK_PATH,
                       admx_exec  = DEFAULT_ADMIXTURE_PATH) {

    message("ADMIXTURE PIPELINE — started: ", Sys.time())
    if (!file.exists(input_vcf)) stop("Input VCF not found: ", input_vcf, call. = FALSE)
    validate_sample_info(sample_info)
    check_dir_exists_or_create(output_dir)

    input_plink   <- make_bed(input_vcf, output_dir, sample_info, exec_plink)
    pop_files     <- make_pop_file(input_plink, sample_info)
    pred_admixture <- run_ADMIXTURE(input_plink, pop_files, output_dir, admx_exec)

    message("ADMIXTURE PIPELINE — completed: ", Sys.time())
    if (!is.na(pred_admixture$output_pop))
        message("  Population level    : ", pred_admixture$output_pop)
    if (!is.na(pred_admixture$output_subpop))
        message("  Sub-population level: ", pred_admixture$output_subpop)
    return(pred_admixture)
}


## ============================================================
## Main Entry Point: ADMIXTURE with Optional SNP Filtering
## ============================================================

#' Run ADMIXTURE with optional discriminant SNP pre-filtering
#'
#' Main entry point called by admixed_inference.r. If a raw VCF is
#' provided (with 'CHR' placeholder), first runs extract_discriminant_SNPs()
#' to select the n_snps most informative markers (ranked by a pre-computed
#' ancestry discrimination score), then runs the full ADMIXTURE pipeline.
#' If a filtered_vcf is provided directly, the SNP filtering step is skipped.
#'
#' @param samples_info  Character. Path to sample info TSV.
#' @param output_dir    Character. Root output directory.
#' @param raw_vcf       Character or NA. Per-chromosome VCF template with
#'                      'CHR' placeholder. Supply this OR filtered_vcf.
#' @param filtered_vcf  Character or NA. Path to pre-filtered single VCF.
#' @param n_snps        Integer or NA. Number of discriminant SNPs to retain
#'                      (required if raw_vcf supplied).
#' @param maf           Numeric or NA. Minor allele frequency threshold
#'                      (required if raw_vcf supplied).
#' @param score_file    Character or NA. Path to per-SNP discrimination scores
#'                      (required if raw_vcf supplied).
#' @param threads       Integer. Parallel threads for SNP extraction (default 10).
#' @param subpop        Character or NULL. If non-null, restrict ancestry
#'                      filtering to this sub-population set.
#' @return Named list with output_pop and output_subpop prediction file paths.
get_admixture_prediction <- function(samples_info, output_dir,
                                      raw_vcf = NA, filtered_vcf = NA,
                                      n_snps = NA, maf = NA, score_file = NA,
                                      threads = 10, subpop = NULL) {

    ## -- Validate sample info --
    check_file_exists(samples_info, "samples_info")
    samples_info_df <- read.table(samples_info, header = TRUE, sep = "\t")
    validate_sample_info(samples_info_df)
    check_dir_exists_or_create(output_dir)

    ## -- Validate VCF input: exactly one of raw_vcf or filtered_vcf --
    if (!is.na(raw_vcf) && !is.na(filtered_vcf))
        stop("Provide only one VCF input: raw_vcf OR filtered_vcf", call. = FALSE)
    if (is.na(raw_vcf) && is.na(filtered_vcf))
        stop("Provide at least one VCF input: raw_vcf or filtered_vcf", call. = FALSE)

    input_vcf <- NULL

    if (!is.na(raw_vcf)) {
        ## raw_vcf must use 'CHR' as chromosome placeholder
        if (!grepl("CHR", raw_vcf))
            stop("raw_vcf must include 'CHR' as chromosome placeholder", call. = FALSE)
        chr_paths  <- sapply(1:22, function(chr) gsub("CHR", chr, raw_vcf))
        chr_exists <- file.exists(chr_paths)
        input_vcf  <- chr_paths[chr_exists]
        missing    <- which(!chr_exists)
        if (length(input_vcf) == 0) stop("No VCF files found for chr 1–22", call. = FALSE)
        if (length(missing) > 0)    warning("Missing VCF for chromosomes: ", paste(missing, collapse = ", "))

        ## SNP filtering parameters are mandatory when starting from raw VCF
        if (is.na(n_snps))     stop("n_snps is required with raw_vcf",    call. = FALSE)
        if (is.na(maf))        stop("maf is required with raw_vcf",        call. = FALSE)
        if (is.na(score_file)) stop("score_file is required with raw_vcf", call. = FALSE)
        check_file_exists(score_file, "score_file")
    }

    if (!is.na(filtered_vcf)) {
        check_file_exists(filtered_vcf, "filtered_vcf")
        input_vcf <- filtered_vcf
    }

    ## -- Optional Step 0: Extract discriminant SNPs --
    if (!is.na(raw_vcf)) {
        message(">>> Extracting discriminant SNPs...")
        ## Ancestry list from pop column (training samples only), or user-supplied subpop
        ancestries <- if (is.null(subpop)) {
            anc <- na.omit(unique(samples_info_df$pop))
            anc[!anc %in% c("-", "", " ")]
        } else {
            c(subpop)
        }
        input_vcf <- extract_discriminant_SNPs(input_vcf, score_file, samples_info_df,
                                                ancestries, maf, n_snps, output_dir)
        message("  -> Discriminant SNP extraction complete.")
    }

    ## -- Steps 1–3: ADMIXTURE --
    message(">>> Running supervised ADMIXTURE...")
    output_admixture <- ADMIXTURE(input_vcf, samples_info_df, output_dir)
    message(">>> get_admixture_prediction complete.")
    return(output_admixture)
}