# ---- Helper: décide quel haplotype porte le SNV ----
which_phase <- function(phase, SNV_with_SNP) {
  if (is.na(phase) || is.na(SNV_with_SNP)) return(NA)
  if ((phase == "0|1" &&  SNV_with_SNP) || (phase == "1|0" && !SNV_with_SNP)) {
    return("hap1")
  } else if ((phase == "1|0" &&  SNV_with_SNP) || (phase == "0|1" && !SNV_with_SNP)) {
    return("hap2")
  } else {
    return(NA)
  }
}
check_hz_haplotype <- function(pred_SNV, pred_hap1, pred_hap2, hap_SNV, hapA, hapB) {
  
  if (length(pred_SNV) != 2 || length(unique(pred_SNV)) != 2) {
    return(c(ancestry_REF = NA, ancestry_SNV = NA, haplotype_SNV = NA))
  }
  
  # Cas 1 : mapping direct
  if (pred_SNV[1] %in% hapA && pred_SNV[2] %in% hapB) {
    hap_map <- c(hap1 = "hapA", hap2 = "hapB")
    
  # Cas 2 : mapping inversé
  } else if (pred_SNV[2] %in% hapA && pred_SNV[1] %in% hapB) {
    hap_map <- c(hap1 = "hapB", hap2 = "hapA")
    
  } else {
    return(c(ancestry_REF = NA, ancestry_SNV = NA, haplotype_SNV = NA))
  }
  
  ancestry_SNV   <- pred_SNV[[hap_SNV]]
  haplotype_SNV  <- hap_map[[hap_SNV]]
  ancestry_REF = setdiff(c(pred_hap1,pred_hap2), ancestry_SNV)
  c(ancestry_REF = ancestry_REF, ancestry_SNV = ancestry_SNV, haplotype_SNV = haplotype_SNV)
}

# ---- Inférence d'ancestry pour un chromosome ----
get_ancestry_chr <- function(sample, hapA, hapB, chr, pred_gnomix, phased_SNV_chr,
                             input_file_beagle, pop_vector, verbose = FALSE) {

    cols <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
            "INFO", "FORMAT", sample)
    SNP_sample <- fread(input_file_beagle, skip = "#CHROM", sep = "\t",
                      header = TRUE, select = cols, data.table = FALSE)

    hap1_col <- paste0(sample, ".0")
    hap2_col <- paste0(sample, ".1")
    pred_gnomix$pos_seg <- (pred_gnomix$spos + pred_gnomix$epos) / 2

    prediction <- lapply(seq_len(nrow(phased_SNV_chr)), function(idx) {
        row <- as.vector(unlist(phased_SNV_chr[idx, ]))
        chrs          <- as.integer(unlist(strsplit(row[4], split = ",")))
        pos_snps      <- as.integer(unlist(strsplit(row[5], split = ",")))
        phasing_state <- as.logical(unlist(strsplit(row[length(row)], split = ",")))

        genotype <- SNP_sample[SNP_sample$POS %in% pos_snps, sample]
        if (length(genotype) == 0) return(c(ancestry_REF = NA, ancestry_SNV = NA, haplotype_SNV = NA))

        # --- Filtre accuracy : on ne garde que les SNVs phasés à 100% ---
        cont_table <- table(phasing_state, genotype)
        if (length(cont_table) == 0 || sum(cont_table) == 0) return(c(ancestry_REF = NA, ancestry_SNV = NA, haplotype_SNV = NA))
        main_diag_sum      <- sum(diag(cont_table))
        secondary_diag_sum <- sum(diag(cont_table[, ncol(cont_table):1, drop = FALSE]))
        accuracy <- max(main_diag_sum / sum(cont_table),
                        secondary_diag_sum / sum(cont_table))
        if (accuracy != 1) return(c(ancestry_REF = NA, ancestry_SNV = NA, haplotype_SNV = NA))

        # Position représentative
        if (length(pos_snps) > 1) {
          pos_snp <- pos_snps[which.min(abs(pos_snps - median(pos_snps)))]
        } else {
          pos_snp <- pos_snps
        }

        hap_SNV <- which_phase(genotype[1], phasing_state[1])
        if (is.na(hap_SNV)) return(c(ancestry_REF = NA, ancestry_SNV = NA, haplotype_SNV = NA))

        # Segment gnomix le plus proche
        closest <- which.min(abs(pred_gnomix$pos_seg - pos_snp))
        pred_hap <- as.vector(unlist(pred_gnomix[closest, c(hap1_col, hap2_col)]))

        pred_hap1 <- pop_vector[pred_hap[1] + 1]   # +1 si codes commencent à 0
        pred_hap2 <- pop_vector[pred_hap[2] + 1]
        pred_SNV <- c(hap1 = pred_hap1, hap2 = pred_hap2)
        return(check_hz_haplotype(pred_SNV, pred_hap1, pred_hap2, hap_SNV, hapA, hapB))
    })


    prediction_df <- as.data.frame( do.call(rbind, prediction), stringsAsFactors = FALSE)
        
    phased_SNV_chr = cbind(phased_SNV_chr, prediction_df)
    return(phased_SNV_chr)
}

# ---- Inférence d'ancestry pour un sample (tous chr) ----
get_ancestry <- function(sample, hapA, hapB, chrs = 1:22, phased_SNV_path, input_pred,
                         input_beagle, output_path, verbose = FALSE, rewrite=FALSE) {

  message(">>> SNVs ancestry inference for sample: ", sample)

  if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)
    
  ancestry_SNV_path <- file.path(output_path, paste0("SNV_ancestry_", sample, ".txt"))

  if (file.exists(ancestry_SNV_path) & !rewrite) {
    message("   -> Already exists: ", ancestry_SNV_path)
    return(ancestry_SNV_path)
  }

  if (!file.exists(phased_SNV_path)) {
    warning("Missing phased SNV file for ", sample, ": ", phased_SNV_path)
    return(NULL)
  }
  if (!dir.exists(input_pred)) stop("Missing gnomix directory: ", input_pred)

  missing_beagle <- unlist(lapply(chrs, function(chr) {
    f <- sub("chrCHR", paste0("chr", chr), input_beagle)
    if (!file.exists(f)) f else NULL
  }))
  if (length(missing_beagle) > 0) {
    warning("Missing Beagle files for ", sample, ":\n",
            paste(missing_beagle, collapse = "\n"))
    return(NULL)
  }

  phased_SNV <- fread(phased_SNV_path, header = TRUE, sep = "\t",
                      data.table = FALSE)

  pred_SNV_list <- mclapply(chrs, function(chr) {
    pred_file <- file.path(input_pred, paste0("chr", chr), "query_results.msp")
    if (!file.exists(pred_file)) {
      warning("Missing gnomix prediction: ", pred_file)
      return(NULL)
    }

    line1 <- readLines(pred_file, n = 1)
    line1 <- sub("^#Subpopulation order/codes:\\s*", "", line1)
    pops  <- strsplit(line1, "\t")[[1]]
    pop_vector <- sub("=.*", "", pops)

    pred_gnomix <- fread(pred_file, header = TRUE, data.table = FALSE)
    phased_SNV_chr <- phased_SNV[phased_SNV$chr_snv == chr, ]

    if (nrow(phased_SNV_chr) == 0) return(NULL)

    input_file_beagle <- gsub("chrCHR", paste0("chr", chr), input_beagle)
    get_ancestry_chr(sample,  hapA, hapB, chr, pred_gnomix, phased_SNV_chr, input_file_beagle, pop_vector, verbose)

  },mc.cores=11)

  pred_SNV_list <- Filter(Negate(is.null), pred_SNV_list)
  if (length(pred_SNV_list) == 0) {
    warning("No predictions produced for ", sample)
    return(NULL)
  }

  pred_SNV_sample <- as.data.frame(do.call(rbind, pred_SNV_list))
  pred_SNV_sample <- pred_SNV_sample[order(pred_SNV_sample$chr_snv,
                                           pred_SNV_sample$pos_snv), ]

  write.table(pred_SNV_sample, file = ancestry_SNV_path,
              quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  message(">>> Done: ", sample)
  return(ancestry_SNV_path)
}