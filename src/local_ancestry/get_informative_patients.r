library(data.table)
library(igraph)
library(GenomicRanges)
library(karyoploteR)
library(RColorBrewer)

SAMPLE_INFO_PATH <- "/srv/home/mlef0011/Phasomix/rawdata/1kGp_PCAWG/info/samples_info_no_AMR_1kGp_train_PCAWG_test.txt"
samples_info <- read.table(SAMPLE_INFO_PATH, header = TRUE, sep = "\t")



# ============================================================
# Fonction : genere un karyoplot pour un patient
# ============================================================
# Arguments :
#   sid              : sample ID
#   pred_all_chr     : df des predictions Gnomix (toutes les fenetres)
#   ancestry         : nom du niveau d'ancestry (pour le titre)
#   cols_ancestries  : vecteur nomme de couleurs par ancestry
#   cna_base_icgc    : chemin vers le dossier CNA ICGC
#   cna_base_tcga    : chemin vers le dossier CNA TCGA
#
# Retourne : invisible(NULL), produit un plot sur le device graphique actif
# ============================================================
get_karyoplot <- function(sid,
                          pred_all_chr,
                          ancestry,
                          cols_ancestries,
                          cna_base_icgc = "/srv/home/mlef0011/rawdata/PCAWG_data/ICGC/cna/",
                          cna_base_tcga = "/srv/home/mlef0011/rawdata/PCAWG_data/TCGA/cna/") {

   # ---- 1. localiser le fichier CNA ----
   cna_path1 <- paste0(cna_base_icgc, sid,
                       ".consensus.20170119.somatic.cna.annotated.txt")
   cna_path2 <- paste0(cna_base_tcga, sid,
                       ".consensus.20170119.somatic.cna.txt")
   cna_path <- if (file.exists(cna_path1)) {
      cna_path1
   } else if (file.exists(cna_path2)) {
      cna_path2
   } else {
      NA_character_
   }

   if (is.na(cna_path)) {
      warning("No CNA file for sample ", sid)
      return(invisible(NULL))
   }

   # ---- 2. lire et filtrer les zones d'imbalance ----
   cna <- read.table(cna_path, header = TRUE, sep = "\t")
   iz <- cna[cna$chromosome %in% 1:22 &
                !is.na(cna$major_cn) &
                !is.na(cna$minor_cn) &
                cna$major_cn != cna$minor_cn, ]

   if (nrow(iz) > 0) {
      iz$chromosome <- paste0("chr", iz$chromosome)
      iz <- iz[order(iz$chromosome, iz$start), ]
      iz$length <- iz$end - iz$start + 1
      iz <- iz[iz$length > 1e7, ]
      iz$color <- rep(c("black", "yellow"), length.out = nrow(iz))
      gr_iz <- GRanges(
         seqnames = iz$chromosome,
         ranges   = IRanges(start = iz$start, end = iz$end),
         color    = iz$color
      )
   } else {
      gr_iz <- GRanges()
   }

   # ---- 3. preparer le df ancestry pour ce patient ----
   df <- pred_all_chr[, c(colnames(pred_all_chr)[1:6],
                          paste0(sid, c(".0", ".1")))]
   colnames(df)[(ncol(df)-1):ncol(df)] <- c("pred_phase1", "pred_phase2")
   df$`#chm` <- paste0("chr", df$`#chm`)

   tbl <- table(c(df$pred_phase1, df$pred_phase2))
   pct <- round(100 * tbl / sum(tbl), 2)

   gr <- GRanges(seqnames    = df$`#chm`,
                 ranges      = IRanges(start = df$spos, end = df$epos),
                 pred_phase1 = df$pred_phase1,
                 pred_phase2 = df$pred_phase2)

   # ---- 4. plot ----
   kp <- plotKaryotype(genome      = "hg19",
                       chromosomes = unique(df$`#chm`),
                       cex         = 0.7)

   kpAddMainTitle(kp,
                  paste0("Local ancestry\n", sid, " (", ancestry, ")"),
                  cex = 1)

   # zones d'imbalance
   for (color_val in unique(gr_iz$color)) {
      regions_color <- gr_iz[gr_iz$color == color_val]
      col <- adjustcolor(color_val, alpha.f = 0.15)
      kpPlotRegions(kp,
                    data   = regions_color,
                    col    = col,
                    r0     = -0.85,
                    r1     = 0.40,
                    border = "#555555",
                    lwd    = 0.5)
   }

   # ancestries par chromosome
   for (chr in unique(seqnames(gr))) {
      chr.dp <- sort(keepSeqlevels(gr, value = chr, pruning.mode = "coarse"))

      # Phase 1
      for (anc in names(cols_ancestries)) {
         kpPlotRegions(kp,
                       data   = chr.dp[chr.dp$pred_phase1 == anc],
                       r0     = 0, r1 = 0.2,
                       col    = cols_ancestries[anc],
                       border = NA)
      }

      # Phase 2
      for (anc in names(cols_ancestries)) {
         kpPlotRegions(kp,
                       data   = chr.dp[chr.dp$pred_phase2 == anc],
                       r0     = -0.45, r1 = -0.65,
                       col    = cols_ancestries[anc],
                       border = NA)
      }
   }

   # legende
   legend_labels <- paste0(names(cols_ancestries),
                           " (", pct[names(cols_ancestries)], "%)")
   legend("bottomright",
          legend = legend_labels,
          fill   = cols_ancestries,
          title  = "",
          cex    = 0.8,
          bg     = "white",
          bty    = "n")

   invisible(NULL)
}

compute_stats <- function(df_pred, col0, col1) {

   # proportions par ancestry (tous haps confondus)
   pops_vec <- c(df_pred[[col0]], df_pred[[col1]])
   pct <- round(100 * prop.table(table(pops_vec)), 2)

   # 2eme ancestry la plus frequente
   pct_sorted <- sort(pct, decreasing = TRUE)
   second_anc_name <- if (length(pct_sorted) >= 2) names(pct_sorted)[2] else NA_character_
   second_anc_prop <- if (length(pct_sorted) >= 2) as.numeric(pct_sorted[2]) else 0

   # proportions par combinaison hap0-hap1 (ordre alphabetique)
   comb <- paste(pmin(df_pred[[col0]], df_pred[[col1]]),
                 pmax(df_pred[[col0]], df_pred[[col1]]),
                 sep = "-")
   pct_comb <- round(100 * prop.table(table(comb)), 2)

   # top combi Hz
   is_hz <- vapply(strsplit(names(pct_comb), "-"),
                   function(x) x[1] != x[2],
                   logical(1))
   if (any(is_hz)) {
      pct_hz <- pct_comb[is_hz]
      top_hz_name <- names(pct_hz)[which.max(pct_hz)]
      top_hz_prop <- as.numeric(max(pct_hz))
   } else {
      top_hz_name <- NA_character_
      top_hz_prop <- 0
   }

   # min_hom_prop
   if (!is.na(top_hz_name)) {
      ancs <- strsplit(top_hz_name, "-")[[1]]
      hom1 <- paste(ancs[1], ancs[1], sep = "-")
      hom2 <- paste(ancs[2], ancs[2], sep = "-")
      p1 <- if (hom1 %in% names(pct_comb)) as.numeric(pct_comb[hom1]) else 0
      p2 <- if (hom2 %in% names(pct_comb)) as.numeric(pct_comb[hom2]) else 0
      min_hom_prop <- min(p1, p2)
   } else {
      min_hom_prop <- NA_real_
   }

   list(second_anc      = second_anc_name,
        second_anc_prop = second_anc_prop,
        top_hz          = top_hz_name,
        top_hz_prop     = top_hz_prop,
        min_hom_prop    = min_hom_prop)
}

# ============================================================
# Helpers
# ============================================================

get_cna_path <- function(sid) {
  p1 <- sprintf("/srv/home/mlef0011/rawdata/PCAWG_data/ICGC/cna/%s.consensus.20170119.somatic.cna.annotated.txt", sid)
  p2 <- sprintf("/srv/home/mlef0011/rawdata/PCAWG_data/TCGA/cna/%s.consensus.20170119.somatic.cna.txt", sid)
  if (file.exists(p1)) return(p1)
  if (file.exists(p2)) return(p2)
  NA_character_
}

load_iz <- function(sid, min_length = 1e7) {
  cna_path <- get_cna_path(sid)
  if (is.na(cna_path)) return(NULL)
  
  cna <- tryCatch(read.table(cna_path, header = TRUE, sep = "\t"), error = function(e) NULL)
  if (is.null(cna)) return(NULL)
  
  iz <- cna[cna$chromosome %in% 1:22 &
            !is.na(cna$major_cn) & !is.na(cna$minor_cn) &
            cna$major_cn != cna$minor_cn, ]
  iz$length <- iz$end - iz$start + 1
  iz <- iz[iz$length > min_length, ]
  iz <- iz[order(iz$chromosome, iz$start), ]
  if (nrow(iz) == 0) return(NULL)
  iz
}

build_pred_iz <- function(iz, pred_all_chr, col0, col1) {
  gr_iz <- GRanges(as.character(iz$chromosome), IRanges(iz$start, iz$end))
  gr_pred <- GRanges(as.character(pred_all_chr$`#chm`),
                     IRanges(pred_all_chr$spos, pred_all_chr$epos))
  
  hits <- findOverlaps(gr_pred, gr_iz)
  if (length(hits) == 0) return(NULL)
  
  inter <- pintersect(gr_pred[queryHits(hits)], gr_iz[subjectHits(hits)])
  dt <- as.data.table(pred_all_chr[queryHits(hits), c(col0, col1)])
  dt[, w := width(inter)]
  dt[, `:=`(a1 = pmin(get(col0), get(col1)),
            a2 = pmax(get(col0), get(col1)))]
  dt
}

# ============================================================
# Bipartition / classification
# ============================================================

try_bipartition <- function(g) {
  n <- vcount(g)
  colors <- rep(NA_integer_, n)
  conflicts <- integer(0)
  
  for (start in seq_len(n)) {
    if (!is.na(colors[start])) next
    queue <- start
    colors[start] <- 0
    while (length(queue) > 0) {
      v <- queue[1]; queue <- queue[-1]
      for (u in as.integer(neighbors(g, v))) {
        if (is.na(colors[u])) {
          colors[u] <- 1 - colors[v]
          queue <- c(queue, u)
        } else if (colors[u] == colors[v]) {
          conflicts <- c(conflicts, get.edge.ids(g, c(v, u)))
        }
      }
    }
  }
  
  if (length(conflicts) == 0) {
    return(list(is_bipartite = TRUE, colors = colors, ambiguous_vertices = character(0)))
  }
  
  bc <- biconnected_components(g)
  ambiguous <- character(0)
  for (i in seq_along(bc$component_edges)) {
    if (any(as.integer(bc$component_edges[[i]]) %in% conflicts)) {
      ambiguous <- union(ambiguous, V(g)$name[as.integer(bc$components[[i]])])
    }
  }
  list(is_bipartite = FALSE, colors = colors,
       ambiguous_vertices = ambiguous, conflict_edges = conflicts)
}

classify_ancestries <- function(g) {
  ancestries <- V(g)$name
  comps <- igraph::components(g)
  result <- data.frame(ancestry = ancestries, component = comps$membership,
                       haplotype = NA_character_, status = NA_character_,
                       stringsAsFactors = FALSE)
  
  for (cid in unique(comps$membership)) {
    sub_g <- induced_subgraph(g, which(comps$membership == cid))
    bip <- try_bipartition(sub_g)
    idx <- match(V(sub_g)$name, result$ancestry)
    
    if (bip$is_bipartite) {
      result$haplotype[idx] <- ifelse(bip$colors == 0,
                                      sprintf("hap_A_comp%d", cid),
                                      sprintf("hap_B_comp%d", cid))
      result$status[idx] <- "resolved"
    } else {
      is_amb <- V(sub_g)$name %in% bip$ambiguous_vertices
      result$status[idx] <- ifelse(is_amb, "ambiguous", "resolved")
      result$haplotype[idx] <- ifelse(is_amb, NA_character_,
                                      ifelse(bip$colors == 0,
                                             sprintf("hap_A_comp%d", cid),
                                             sprintf("hap_B_comp%d", cid)))
    }
  }
  result
}

# ============================================================
# Étape 1 : sélection des ancestries à vérifier
# ============================================================

get_ancestries_to_check <- function(sid, pred_all_chr, meta_cols,
                                    hetero_min = 5,
                                    hm_threshold = 0.1, max_hm_zones = 2) {
  col0 <- paste0(sid, ".0")
  col1 <- paste0(sid, ".1")
  pred_sample <- pred_all_chr[, c(meta_cols, col0, col1)]
  all_ancestries <- sort(unique(c(pred_sample[[col0]], pred_sample[[col1]])))
  
  mat_raw <- prop.table(table(
    factor(pred_sample[[col0]], levels = all_ancestries),
    factor(pred_sample[[col1]], levels = all_ancestries))) * 100

  mat <- mat_raw + t(mat_raw)
  diag(mat) <- diag(mat_raw)

  homo <- diag(mat)
  hetero <- rowSums(mat) - homo
  candidates_hz <- names(homo)[hetero > hetero_min]
  candidates_hm =   candidates_hz[sapply(candidates_hz, function(x) { homo[x]/ (hetero[x]+homo[x]) < hm_threshold})]
  candidates = intersect(candidates_hz, candidates_hm)
  if (length(candidates) == 0) return(character(0))
  
  iz <- load_iz(sid)
  if (is.null(iz)) return(character(0))
  
  gr_pred <- GRanges(as.character(pred_all_chr$`#chm`),
                     IRanges(pred_all_chr$spos, pred_all_chr$epos))
  
  candidates[sapply(candidates, function(ancestry) {
    states <- sapply(seq_len(nrow(iz)), function(i) {
      gr_iz <- GRanges(iz$chromosome[i], IRanges(iz$start[i], iz$end[i]))
      hits <- findOverlaps(gr_pred, gr_iz)
      if (length(hits) == 0) return(NA_character_)
      seg <- pred_all_chr[queryHits(hits), c(col0, col1)]
      p1 <- mean(seg[[col0]] == ancestry, na.rm = TRUE); if (is.na(p1)) p1 <- 0
      p2 <- mean(seg[[col1]] == ancestry, na.rm = TRUE); if (is.na(p2)) p2 <- 0
      if (min(p1, p2) > hm_threshold) "hm" else "hz"
    })
    sum(states == "hm", na.rm = TRUE) <= max_hm_zones
  })]
}

# ============================================================
# Étape 2 : graphe de co-portage et bipartition
# ============================================================

process_sample <- function(sid, pred_all_chr, ancestries_hz,
                           pair_threshold = 1, hom_threshold = 5) {
  anc_check <- ancestries_hz[[sid]]
  if (is.null(anc_check) || length(anc_check) == 0) return(NULL)
  
  iz <- load_iz(sid)
  if (is.null(iz)) return(NULL)
  
  col0 <- paste0(sid, ".0")
  col1 <- paste0(sid, ".1")
  pred_in_iz <- build_pred_iz(iz, pred_all_chr, col0, col1)
  if (is.null(pred_in_iz)) return(NULL)
  
  pair_weights <- pred_in_iz[, .(bp = sum(w)), by = .(a1, a2)]
  pair_weights[, prop := 100 * bp / sum(bp)]
  pair_weights <- pair_weights[order(-prop)]
  
  het_pairs <- pair_weights[a1 != a2 & prop > pair_threshold &
                            a1 %in% anc_check & a2 %in% anc_check]
  hom_pairs <- pair_weights[a1 == a2 & prop > hom_threshold]
  
  classification <- NULL
  conflict_pairs <- NULL
  ambiguous <- character(0)
  
  if (nrow(het_pairs) > 0) {
    g <- graph_from_data_frame(het_pairs[, .(from = a1, to = a2, weight = prop)],
                               directed = FALSE)
    classification <- classify_ancestries(g)
    bip <- try_bipartition(g)
    ambiguous <- classification$ancestry[classification$status == "ambiguous"]
    if (length(bip$conflict_edges) > 0) {
      conflict_pairs <- igraph::as_data_frame(g, what = "edges")[bip$conflict_edges, ]
    }
  }
  
  list(sample = sid,
       ancestries_checked = anc_check,
       het_pairs = het_pairs,
       hom_pairs = hom_pairs,
       pred_in_iz = pred_in_iz,
       col0 = col0,
       col1 = col1,
       classification = classification,
       ambiguous = ambiguous,
       conflict_pairs = conflict_pairs)
}

# ============================================================
# Étape 3 : résumé par patient
# ============================================================

summarize_haplotypes <- function(r) {
  shared <- if (!is.null(r$hom_pairs) && nrow(r$hom_pairs) > 0) {
    as.character(r$hom_pairs$a1)
  } else character(0)
  
  resolved_a <- character(0)
  resolved_b <- character(0)
  if (!is.null(r$classification)) {
    cl <- r$classification[r$classification$status == "resolved", ]
    resolved_a <- cl$ancestry[grepl("^hap_A", cl$haplotype)]
    resolved_b <- cl$ancestry[grepl("^hap_B", cl$haplotype)]
  }
  
  unresolved <- setdiff(r$ancestries_checked,
                        c(shared, resolved_a, resolved_b, r$ambiguous))
  
  unknown <- character(0)
  anchored_a <- character(0)
  anchored_b <- character(0)
  
  if (length(unresolved) > 0) {
    if (length(shared) == 0) {
      unknown <- unresolved
    } else {
      anchor <- shared[1]
      for (anc in unresolved) {
        case1 <- r$pred_in_iz[get(r$col0) == anchor & get(r$col1) == anc, sum(w)]
        case2 <- r$pred_in_iz[get(r$col0) == anc & get(r$col1) == anchor, sum(w)]
        
        if (case1 + case2 == 0) {
          unknown <- c(unknown, anc)
        } else if (case1 > case2) {
          anchored_b <- c(anchored_b, anc)
        } else {
          anchored_a <- c(anchored_a, anc)
        }
      }
    }
  }
  
  hap_A <- sort(unique(c(shared, resolved_a, anchored_a)))
  hap_B <- sort(unique(c(shared, resolved_b, anchored_b)))
  
  status <- if (length(unknown) > 0 || length(r$ambiguous) > 0) "partial" else "resolved"
  
  data.table(
    sample = r$sample,
    status = status,
    hap_A = sprintf("{%s}", paste(hap_A, collapse = ",")),
    hap_B = sprintf("{%s}", paste(hap_B, collapse = ",")),
    n_unknown = length(unknown),
    unknown = paste(unknown, collapse = ","),
    n_ambiguous = length(r$ambiguous),
    ambiguous = paste(r$ambiguous, collapse = ",")
  )
}
                  
chrs <- 1:22
pred_path <- "/srv/home/mlef0011/Phasomix/output/local_ancestry/prediction/high_coverage/IS_corrected/"
ancestry_levels = "POP_no_AMR"

# ============================================================
# Boucle sur chaque ancestry level
# ============================================================
for (ancestry in ancestry_levels) {

   message("\n========== Processing ancestry level: ", ancestry, " ==========\n")

   # ---- Charger pred_all_chr UNE FOIS ----
   all_chr_data <- lapply(chrs, function(chr) {
      pred_file <- paste0(pred_path, ancestry,
                          "/prediction/0.1cM/chr", chr,
                          "/query_results.msp")
      if (file.exists(pred_file)) {
         line1 <- readLines(pred_file, n = 1)
         line1 <- sub("^#Subpopulation order/codes:\\s*", "", line1)
         pops <- strsplit(line1, "\t")[[1]]
         pop_vector <- sub("=.*", "", pops)
         pred_gnomix <- fread(pred_file, header = TRUE, sep = "\t",
                              data.table = FALSE)
         pred_gnomix[, 7:ncol(pred_gnomix)] <- lapply(
            pred_gnomix[, 7:ncol(pred_gnomix)],
            function(col) pop_vector[col + 1]
         )
         return(pred_gnomix)
      }
   })

   all_chr_data <- Filter(Negate(is.null), all_chr_data)
   if (length(all_chr_data) == 0) {
      warning("No prediction files found for ancestry ", ancestry, ", skipping")
      next
   }

   common_cols  <- Reduce(intersect, lapply(all_chr_data, colnames))
   pred_all_chr <- as.data.frame(do.call(rbind, lapply(all_chr_data,
                                                       function(df) df[, common_cols])))
   if (nrow(pred_all_chr) < 1) stop("No prediction data loaded")

   meta_cols <- c("#chm", "spos", "epos", "sgpos", "egpos", "n snps")
   hap_cols  <- setdiff(colnames(pred_all_chr), meta_cols)
   samples   <- unique(sub("\\.[01]$", "", hap_cols))

   message("  -> ", length(samples), " samples loaded")

   # ============================================================
   # Partie 1 : calcul des stats par patient (results_df)
   # ============================================================
   results_df <- do.call(rbind, lapply(samples, function(sid) {

      col0 <- paste0(sid, ".0")
      col1 <- paste0(sid, ".1")
      if (!all(c(col0, col1) %in% colnames(pred_all_chr))) return(NULL)

      # stats genome-wide
      stats_all <- compute_stats(pred_all_chr, col0, col1)

      # CNA
      cna_path1 <- paste0("/srv/home/mlef0011/rawdata/PCAWG_data/ICGC/cna/",
                          sid, ".consensus.20170119.somatic.cna.annotated.txt")
      cna_path2 <- paste0("/srv/home/mlef0011/rawdata/PCAWG_data/TCGA/cna/",
                          sid, ".consensus.20170119.somatic.cna.txt")
      cna_path <- if (file.exists(cna_path1)) cna_path1
                  else if (file.exists(cna_path2)) cna_path2
                  else NA_character_

      if (is.na(cna_path)) return(NULL)

      cna <- tryCatch(
         read.table(cna_path, header = TRUE, sep = "\t"),
         error = function(e) NULL
      )
      if (is.null(cna)) return(NULL)

      iz <- cna[cna$chromosome %in% 1:22 &
                   !is.na(cna$major_cn) &
                   !is.na(cna$minor_cn) &
                   cna$major_cn != cna$minor_cn, ]
      if (nrow(iz) == 0) return(NULL)

      iz$length <- iz$end - iz$start + 1
      iz <- iz[iz$length > 1e7, ]
      if (nrow(iz) <= 10) return(NULL)

      iz <- iz[order(iz$chromosome, iz$start), ]

      gr_iz   <- GRanges(seqnames = iz$chromosome,
                         ranges   = IRanges(start = iz$start, end = iz$end))
      gr_pred <- GRanges(seqnames = pred_all_chr$`#chm`,
                         ranges   = IRanges(start = pred_all_chr$spos,
                                            end   = pred_all_chr$epos))

      hits <- findOverlaps(gr_pred, gr_iz)
      if (length(hits) == 0) return(NULL)

      pred_in_iz <- pred_all_chr[queryHits(hits), ]

      # stats dans les IZ
      stats_iz <- compute_stats(pred_in_iz, col0, col1)

      # prop globale d'homozygotie dans les IZ
      prop_hom_iz <- round(100 * mean(pred_in_iz[[col0]] == pred_in_iz[[col1]]), 2)

      data.frame(sample             = sid,
                 second_anc         = stats_all$second_anc,
                 second_anc_prop    = stats_all$second_anc_prop,
                 top_hz             = stats_all$top_hz,
                 top_hz_prop        = stats_all$top_hz_prop,
                 min_hom_prop       = stats_all$min_hom_prop,
                 second_anc_iz      = stats_iz$second_anc,
                 second_anc_prop_iz = stats_iz$second_anc_prop,
                 top_hz_iz          = stats_iz$top_hz,
                 top_hz_prop_iz     = stats_iz$top_hz_prop,
                 min_hom_prop_iz    = stats_iz$min_hom_prop,
                 prop_hom_iz        = prop_hom_iz,
                 stringsAsFactors = FALSE)
   }))

   # sauvegarder le results_df
   write.table(results_df,
               file = paste0(pred_path, "stats_local_ancestry_", ancestry, ".tsv"),
               sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

   message("  -> Stats saved (", nrow(results_df), " samples)")

   # ============================================================
   # Partie 2 : filtres case1 / case3 et karyoplots
   # ============================================================
    true_case1 <- results_df[results_df$min_hom_prop_iz < 1 &
                                   results_df$top_hz_prop_iz >= 90, ]

    true_case3 <- results_df[results_df$top_hz_prop_iz > 20 & 
                             results_df$top_hz_prop_iz < 90 & 
                             results_df$min_hom_prop_iz < 1, ]

   F1_like        <- true_case1$sample
   mono_x_admixed <- true_case3$sample

   message("  -> F1_like (case 1): ", length(F1_like))
   message("  -> mono_x_admixed (case 3): ", length(mono_x_admixed))

   # ---- palette auto-generee pour ce level ----
   pred_file <- paste0(pred_path, ancestry, "/prediction/0.1cM/chr22/query_results.msp")
   line1 <- readLines(pred_file, n = 1)
   line1 <- sub("^#Subpopulation order/codes:\\s*", "", line1)
   pops <- strsplit(line1, "\t")[[1]]
   pops_present <- sub("=.*", "", pops)

   if (length(pops_present) <= 9) {
      cols_ancestries <- setNames(brewer.pal(max(3, length(pops_present)), "Set1"),
                                   pops_present)
   } else {
      cols_ancestries <- setNames(
         colorRampPalette(brewer.pal(9, "Set1"))(length(pops_present)),
         pops_present
      )
   }

   # ---- karyoplots F1_like ----
   if (length(F1_like) > 0) {
      output_plot <- paste0(pred_path, "F1_like_PCAWG_karyoplot_", ancestry, ".pdf")
      pdf(output_plot, width = 8, height = 11)
      for (sid in F1_like) {
         print(sid)
         get_karyoplot(sid             = sid,
                       pred_all_chr    = pred_all_chr,
                       ancestry        = ancestry,
                       cols_ancestries = cols_ancestries)
      }
      dev.off()
      message("  -> F1_like karyoplots saved")
   }

   # ---- karyoplots mono_x_admixed ----
   if (length(mono_x_admixed) > 0) {
      output_plot <- paste0(pred_path, "mono_x_admixed_PCAWG_karyoplot_", ancestry, ".pdf")
      pdf(output_plot, width = 8, height = 11)
      for (sid in mono_x_admixed) {
         print(sid)
         get_karyoplot(sid             = sid,
                       pred_all_chr    = pred_all_chr,
                       ancestry        = ancestry,
                       cols_ancestries = cols_ancestries)
      }
      dev.off()
      message("  -> mono_x_admixed karyoplots saved")
   }
                                                       

    # ============================================================
    # Pipeline
    # ============================================================
    admixed = unique(c(F1_like, mono_x_admixed))
    ancestries_hz <- lapply(admixed, get_ancestries_to_check,
                            pred_all_chr = pred_all_chr, meta_cols = meta_cols)
    names(ancestries_hz) <- admixed
    results <- lapply(admixed, process_sample,
                      pred_all_chr = pred_all_chr,
                      ancestries_hz = ancestries_hz)
    names(results) <- admixed
    results <- results[!sapply(results, is.null)]
    summary_df <- rbindlist(lapply(results, summarize_haplotypes))
    write.table(summary_df, paste0(pred_path,"info_F1_like_mono_x_admixed_gnomix_PCAWG_",ancestry,".tsv"), col.names=T,row.names=F, sep = '\t', quote=F)
}

message("\n========== Done. ==========\n")