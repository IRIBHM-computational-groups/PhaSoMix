suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(copynumber)
 
})
 
doPCF <- function(baf, chr, minEvents, gamma)
{
  pos=1:length(baf)
  res <- pcf(data.frame(chr=rep(chr,length(baf)),
                        positions=pos,
                        sample1=baf),
             verbose=F,
             kmin=minEvents,
            gamma=gamma)
  res
}  
                                  
                                                    
compute_SE_IZ <- function(segment_CN, 
                                  segment_beagle,
                                  chr,
                                  minor_cn, major_cn, total_cn,
                                  gamma,
                                  minEvents,
                                  min_run) {
  
    # ---------- STEP 1: PCF-based block flipping ----------
    pcf <- doPCF(segment_CN$BAF_phase1, chr = chr, minEvents = minEvents, gamma = gamma)
    pcf_switch <- subset(pcf, mean < 0.5)
 
    in_pcf_switch <- rep(FALSE, nrow(segment_CN))
    if (nrow(pcf_switch) > 0) {
        for (r in seq_len(nrow(pcf_switch))) {
            start_pos <- pcf_switch$start.pos[r]
            end_pos   <- pcf_switch$end.pos[r]
            in_pcf_switch[start_pos:end_pos] <- TRUE
        }
    }
 
    # Initialisation des colonnes bloc
    segment_CN$phase1_corrected <- segment_CN$phase1
    segment_CN$phase2_corrected <- segment_CN$phase2
    segment_CN$BAF_phase1_corrected <- segment_CN$BAF_phase1
    segment_CN$BAF_phase2_corrected <- segment_CN$BAF_phase2
  
    if (any(in_pcf_switch)) {
        rle_mask <- rle(in_pcf_switch)
        ends <- cumsum(rle_mask$lengths)
        starts <- c(1, head(ends, -1) + 1)
        for (i in seq_along(rle_mask$values)) {
          if (rle_mask$values[i]) {
                idx_correction <- starts[i]:ends[i]
                # swap block
                tmp1 <- segment_CN$phase1_corrected[idx_correction]
                segment_CN$phase1_corrected[idx_correction] <- segment_CN$phase2_corrected[idx_correction]
                segment_CN$phase2_corrected[idx_correction] <- tmp1
                tmpb1 <- segment_CN$BAF_phase1_corrected[idx_correction]
                segment_CN$BAF_phase1_corrected[idx_correction] <- segment_CN$BAF_phase2_corrected[idx_correction]
                segment_CN$BAF_phase2_corrected[idx_correction] <- tmpb1
            }
        }
    }
  
    # ---------- STEP 2: SNP-level refinement using LLR ----------
    count1 <- round(segment_CN$BAF_phase1_corrected * segment_CN$Good_depth)
    count2 <- round(segment_CN$BAF_phase2_corrected * segment_CN$Good_depth)
    n_all  <- count1 + count2
 
    llr_vec <- mapply(function(k, n, a1, a2, minor_cn, major_cn, total_cn) {
        p1 <- (major_cn * a1 + minor_cn * a2) / total_cn
        p2 <- (minor_cn * a1 + major_cn * a2) / total_cn
        p1 <- pmin(pmax(p1, 1e-8), 1 - 1e-8)
        p2 <- pmin(pmax(p2, 1e-8), 1 - 1e-8)
        log(p1 / p2)
    }, k = count1, n = n_all, a1 = segment_CN$BAF_phase1_corrected, 
    a2 = segment_CN$BAF_phase2_corrected,
    minor_cn = minor_cn, major_cn = major_cn, total_cn = total_cn)
 
    segment_CN$llr_smooth <- llr_vec
 
    # SNPs Ã  inverser
    candidate_snps <- which(segment_CN$llr_smooth < 0)
    
    if (length(candidate_snps) > 0){
        cand_rle <- rle(seq_along(segment_CN$llr_smooth) %in% candidate_snps)
        ends <- cumsum(cand_rle$lengths)
        starts <- c(1, head(ends, -1) + 1)
        for (i in seq_along(cand_rle$values)) {
            if (cand_rle$values[i]) {
                idx_run <- starts[i]:ends[i]
                if (length(idx_run) >= min_run) {
                    tmp1 <- segment_CN$phase1_corrected[idx_run]
                    segment_CN$phase1_corrected[idx_run] <- segment_CN$phase2_corrected[idx_run]
                    segment_CN$phase2_corrected[idx_run] <- tmp1
                    tmpb1 <- segment_CN$BAF_phase1_corrected[idx_run]
                    segment_CN$BAF_phase1_corrected[idx_run] <- segment_CN$BAF_phase2_corrected[idx_run]
                    segment_CN$BAF_phase2_corrected[idx_run] <- tmpb1
                }
            }
        }
    }    

    segment_CN$switch <- segment_CN$phase1 != segment_CN$phase1_corrected
 
    rle_mask <- rle(segment_CN$switch)
    ends <- cumsum(rle_mask$lengths)
    starts <- c(1, head(ends, -1) + 1)
    switch_pos <- segment_CN$POS[segment_CN$switch == TRUE]
    
    idx <- segment_beagle$POS %in% switch_pos
    segment_beagle[idx, 10] <- ifelse(
        segment_beagle[idx, 10] == "0|1", "1|0",
        ifelse(segment_beagle[idx, 10] == "1|0", "0|1", segment_beagle[idx, 10])
    )
 
 
    return(list(candidata_snps = candidate_snps,
                segment_CN = segment_CN,
                segment_beagle_corrected = segment_beagle,
                pos_to_switch = switch_pos,
                pos_SE = segment_CN$POS[starts[-1]],
                rle_mask = rle_mask$lengths))
 
 
}