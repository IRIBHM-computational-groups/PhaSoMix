#!/usr/bin/env Rscript

# ============================================================
# Title:     Utility Functions for Ancestry Analysis
# Author:    Maxime Lefebvre
# Created:   2025-11-26
# Purpose:   Common utility functions for file handling,
#            ancestry analysis, and data validation
# Repository: https://github.com/yourusername/Phasomix
# ============================================================

suppressPackageStartupMessages({
    library(data.table)     
    library(dplyr)
})
## ============================================================
## File System Utilities
## ============================================================

#' Check if a file exists and stop execution if not
check_file_exists <- function(path, label) {
  if (!file.exists(path)) {
    stop(sprintf("The %s file does not exist: %s", label, path), call. = FALSE)
  }
  invisible(TRUE)
}


#' Create directory if it doesn't exist
check_dir_exists_or_create <- function(path) {
  if (!dir.exists(path)) {
    message(sprintf("Creating output directory: %s", path))
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  invisible(TRUE)
}

check_all_chr <- function(path) {
  if (!grepl("CHR", path)) {
    stop(paste0(path, " must include 'CHR' as chromosome placeholder."))
  }
  chr_paths <- sapply(1:22, function(chr) gsub("CHR", chr, path))
  chr_exists <- file.exists(chr_paths)
  input_vcf <- chr_paths[chr_exists]
  missing_chr <- which(!chr_exists)
  existing_chrs <- which(chr_exists)

  if (length(input_vcf) == 0) {
    stop(paste0("No ", path, " found for chromosomes 1 to 22."))
  }
  if (length(missing_chr) > 0) {
    warning(sprintf("Missing %s for chromosomes: %s", path, paste(missing_chr, collapse = ", ")))
  }
  return(existing_chrs)
}
                      
## ============================================================
## Ancestry Analysis Utilities
## ============================================================

#' Check if an individual is admixed based on ancestry proportions
#' An individual is classified as admixed if:
#'   1. Sum of top 2 ancestry components > 0.8
#'   2. Absolute difference between top 2 components < 0.2
check_conditions <- function(row) {
  if (length(row) < 2) {
    warning("Row has fewer than 2 values. Cannot determine admixture.")
    return(FALSE)
  }
  
  # Get top 2 values
  top2 <- sort(row, decreasing = TRUE)[1:2]
  
  # Check conditions
  sum_condition <- sum(top2) > 0.8
  diff_condition <- abs(diff(top2)) < 0.2
  
  return(sum_condition & diff_condition)
}


#' Get names of top 2 ancestry components
get_top2_colnames <- function(row) {
  if (length(row) < 2) {
    warning("Row has fewer than 2 values. Returning NA.")
    return(NA_character_)
  }
  
  if (is.null(names(row))) {
    warning("Row has no names. Cannot extract ancestry labels.")
    return(NA_character_)
  }
  
  # Get indices of top 2 values
  top2_indices <- order(row, decreasing = TRUE)[1:2]
  
  # Get corresponding names
  top2_names <- names(row)[top2_indices]
  
  # Sort alphabetically for consistency
  sorted_names <- sort(top2_names)
  
  return(paste(sorted_names, collapse = "-"))
}


#' Get name of top 1 ancestry component
get_top1_colname <- function(row) {
  if (length(row) == 0) {
    warning("Empty row. Returning NA.")
    return(NA_character_)
  }
  
  if (is.null(names(row))) {
    warning("Row has no names. Cannot extract ancestry label.")
    return(NA_character_)
  }
  
  # Get index of maximum value
  top1_index <- which.max(row)
  
  return(names(row)[top1_index])
}
                      
#' check if maximum prediction probability value is higher thant the threshold
check_if_pure <- function(row, tsh) {

  max_val <- max(row, na.rm = TRUE)

  if (max_val >= tsh) {
    return(names(row)[which.max(row)])
  } else {
    return(FALSE)
  }
}


#' Generate all valid ancestry combinations
#' Creates all pairwise combinations of populations and subpopulations
#' within populations.

generate_ancestry_combinations <- function(pop_info) {
  
  message(">>> Generating ancestry combinations...")
  
  # Population-level combinations
  pops_sorted <- sort(unique(pop_info$Super.Population))
  pop_combos <- combn(pops_sorted, 2, simplify = FALSE)
  pop_combinations <- sapply(pop_combos, function(x) paste(sort(x), collapse = "-"))
  
  message(sprintf("  -> Population combinations: %d", length(pop_combinations)))
  
  # Subpopulation combinations (within each super-population)
  subpop_combinations <- unlist(lapply(unique(pop_info$Super.Population), function(pop) {
    
    subpops <- pop_info %>%
      filter(Super.Population == pop) %>%
      pull(Population) %>%
      unique() %>%
      sort()
    
    if (length(subpops) < 2) return(NULL)
    
    combos <- combn(subpops, 2, simplify = FALSE)
    combo_names <- sapply(combos, function(x) paste(sort(x), collapse = "-"))
    
    return(combo_names)
  }))
  
  message(sprintf("  -> Subpopulation combinations: %d", length(subpop_combinations)))
  
  all_combinations <- c(pop_combinations, subpop_combinations)
  
  message(sprintf("  -> Total combinations: %d", length(all_combinations)))
  
  return(all_combinations)
}

                                    

## ============================================================
## Data Validation Utilities
## ============================================================

#' Validate sample information data frame
validate_sample_info <- function(samples_info) {
  
  # Must be a data frame
  if (!is.data.frame(samples_info)) {
    stop("samples_info must be a data frame", call. = FALSE)
  }
  
  # Mandatory columns
  mandatory_cols <- c("sample", "category")
  missing_mandatory <- setdiff(mandatory_cols, colnames(samples_info))
  if (length(missing_mandatory) > 0) {
    stop(
      sprintf(
        "samples_info is missing mandatory columns: %s",
        paste(missing_mandatory, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  # Must have at least one of pop or subpop
  has_pop    <- "pop"    %in% colnames(samples_info)
  has_subpop <- "subpop" %in% colnames(samples_info)
  
  if (!has_pop && !has_subpop) {
    stop("samples_info must contain at least one of: 'pop' or 'subpop'", 
         call. = FALSE)
  }
  
  # Check category values
  if (!all(samples_info$category %in% c(0, 1))) {
    stop("category must contain only 0 (train) or 1 (test)", call. = FALSE)
  }
  
  # NA checks only for columns present
  for (col in intersect(colnames(samples_info), c("sample", "category", "pop", "subpop"))) {
    if (any(is.na(samples_info[[col]]))) {
      warning(sprintf("Column '%s' contains NA values", col), call. = FALSE)
    }
  }
  
  invisible(TRUE)
}

                                       
check_required_columns <- function(df, required_cols) {
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) != 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  df <- df[, c(required_cols, setdiff(colnames(df), required_cols))]
  return(df)
}