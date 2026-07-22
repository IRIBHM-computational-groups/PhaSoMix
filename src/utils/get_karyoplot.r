get_karyoplot <- function(sample, label, pred_path, snv_ancestry_path, info_ma, output_dir) {

   # ------------------------------------------------------------
   # 1. Lire les predictions Gnomix (.msp) chr1-22
   # ------------------------------------------------------------
   all_chr_data <- lapply(1:22, function(chr) {
      pred_file <- paste0(pred_path, label,
                          "/prediction/0.1cM/chr", chr,
                          "/query_results.msp")
      if (file.exists(pred_file)) {
         line1      <- readLines(pred_file, n = 1)
         line1      <- sub("^#Subpopulation order/codes:\\s*", "", line1)
         pops       <- strsplit(line1, "\t")[[1]]
         pop_vector <- sub("=.*", "", pops)
         cols       <- c("#chm", "spos", "epos", "sgpos", "egpos", "n snps",
                         paste0(sample, c(".0", ".1")))
         pred_gnomix <- fread(pred_file, select = cols, header = TRUE, sep = "\t",
                              data.table = FALSE)
         pred_gnomix[, 7:ncol(pred_gnomix)] <- lapply(
            pred_gnomix[, 7:ncol(pred_gnomix)],
            function(col) pop_vector[col + 1]
         )
         return(pred_gnomix)
      } else {
         message("MISSING FILE: ", pred_file)
         return(NULL)
      }
   })

   all_chr_data <- Filter(Negate(is.null), all_chr_data)
   if (length(all_chr_data) == 0) {
      warning("No prediction files found for ", sample, " (", label, "), skipping")
      return(invisible(NULL))
   }

   common_cols  <- Reduce(intersect, lapply(all_chr_data, colnames))
   pred_all_chr <- as.data.frame(do.call(rbind, lapply(all_chr_data,
                                                       function(d) d[, common_cols])))

   # ------------------------------------------------------------
   # 2. Palette ancestry (a partir des populations de chr22)
   # ------------------------------------------------------------
   pred_file    <- paste0(pred_path, label, "/prediction/0.1cM/chr22/query_results.msp")
   line1        <- readLines(pred_file, n = 1)
   line1        <- sub("^#Subpopulation order/codes:\\s*", "", line1)
   pops         <- strsplit(line1, "\t")[[1]]
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

   # ------------------------------------------------------------
   # 3. Zones d'imbalance (CNA major_cn != minor_cn, > 10 Mb)
   # ------------------------------------------------------------
   cna_path <- info_ma[info_ma$sample == sample, "input_CNA"]
   cna      <- read.table(cna_path, header = TRUE, sep = "\t")
   iz       <- cna[cna$chromosome %in% 1:22 &
                      !is.na(cna$major_cn) &
                      !is.na(cna$minor_cn) &
                      cna$major_cn != cna$minor_cn, ]

   if (nrow(iz) > 0) {
      iz$chromosome <- paste0("chr", iz$chromosome)
      iz            <- iz[order(iz$chromosome, iz$start), ]
      iz$length     <- iz$end - iz$start + 1
      iz            <- iz[iz$length > 1e7, ]
   }

   if (nrow(iz) > 0) {
      iz$color <- "black"
      gr_iz <- GRanges(
         seqnames = iz$chromosome,
         ranges   = IRanges(start = iz$start, end = iz$end),
         color    = iz$color
      )
   } else {
      gr_iz <- GRanges()
   }

   # ------------------------------------------------------------
   # 4. df ancestry -> GRanges des segments
   # ------------------------------------------------------------
   df <- pred_all_chr[, c(colnames(pred_all_chr)[1:6],
                          paste0(sample, c(".0", ".1")))]
   colnames(df)[(ncol(df) - 1):ncol(df)] <- c("pred_phase1", "pred_phase2")
   df$`#chm` <- paste0("chr", df$`#chm`)

   gr <- GRanges(seqnames    = df$`#chm`,
                 ranges      = IRanges(start = df$spos, end = df$epos),
                 pred_phase1 = df$pred_phase1,
                 pred_phase2 = df$pred_phase2)

   # ------------------------------------------------------------
   # 5. SNVs -> GRanges + attribution de phase
   # ------------------------------------------------------------
   snv_ancestry          <- read.table(snv_ancestry_path, header = TRUE, sep = "\t")
   snv_ancestry$chr_snv  <- paste0("chr", snv_ancestry$chr_snv)
   snv_ancestry          <- snv_ancestry[!is.na(snv_ancestry$ancestry_SNV), ]

   gr_snv <- GRanges(seqnames = snv_ancestry$chr_snv,
                     ranges   = IRanges(start = snv_ancestry$pos_snv,
                                        end   = snv_ancestry$pos_snv),
                     ancestry = snv_ancestry$ancestry_SNV)

   hits    <- findOverlaps(gr_snv, gr)
   seg_idx <- subjectHits(hits)
   snv_idx <- queryHits(hits)

   anc_snv <- gr_snv$ancestry[snv_idx]
   p1      <- gr$pred_phase1[seg_idx]
   p2      <- gr$pred_phase2[seg_idx]

   phase <- ifelse(anc_snv == p1 & anc_snv != p2, "phase1",
            ifelse(anc_snv == p2 & anc_snv != p1, "phase2", NA))

   gr_snv$phase          <- NA_character_
   gr_snv$phase[snv_idx] <- phase

   # ------------------------------------------------------------
   # 6. Plot
   # ------------------------------------------------------------
   pdf(file.path(output_dir, paste0("karyoplot_local_ancestry_", sample, ".pdf")),
       width = 8, height = 10)
   on.exit(dev.off(), add = TRUE)

   pp <- getDefaultPlotParams(plot.type = 1)
   pp$data1height    <- 120
   pp$data1inmargin  <- 5
   pp$ideogramheight <- 10
   pp$leftmargin     <- 0.08
   pp$topmargin      <- 20
   pp$bottommargin   <- 30

   # Track layout : tout dans [0, -1] pour eviter la troncature
   R0_S1 <- -0.05;  R1_S1 <- -0.20   # ticks SNV phase 1 (au-dessus de P1)
   R0_P1 <- -0.20;  R1_P1 <- -0.45   # bande ancestry phase 1
   R0_P2 <- -0.50;  R1_P2 <- -0.75   # bande ancestry phase 2
   R0_S2 <- -0.75;  R1_S2 <- -0.90   # ticks SNV phase 2 (sous P2)
   R0_IZ <- -0.05;  R1_IZ <- -0.90   # IZ : colonne de fond couvrant tout

   kp <- plotKaryotype(genome      = "hg19",
                       chromosomes = unique(as.character(seqnames(gr))),
                       plot.params = pp,
                       cex         = 0.7)

   kpAddMainTitle(kp, paste0("Local ancestry - ", sample), cex = 1)

   kpAddLabels(kp, labels = "p1", r0 = R0_P1, r1 = R1_P1, cex = 0.5, label.margin = 0.01)
   kpAddLabels(kp, labels = "p2", r0 = R0_P2, r1 = R1_P2, cex = 0.5, label.margin = 0.01)

   # (a) Zones d'imbalance - fond translucide
   for (color_val in unique(gr_iz$color)) {
      regions_color <- gr_iz[gr_iz$color == color_val]
      kpRect(kp, data = regions_color,
             y0 = 0, y1 = 1, r0 = R0_IZ, r1 = R1_IZ,
             col = adjustcolor(color_val, alpha.f = 0.12), border = NA)
   }

   # (b) Bandes ancestry - un appel vectorise par phase
   for (chr in unique(as.character(seqnames(gr)))) {
      chr.dp <- sort(keepSeqlevels(gr, value = chr, pruning.mode = "coarse"))

      kpRect(kp, data = chr.dp, y0 = 0, y1 = 1, r0 = R0_P1, r1 = R1_P1,
             col = cols_ancestries[chr.dp$pred_phase1], border = NA)
      kpRect(kp, data = chr.dp, y0 = 0, y1 = 1, r0 = R0_P2, r1 = R1_P2,
             col = cols_ancestries[chr.dp$pred_phase2], border = NA)
   }

   # (c) Ticks SNV
   for (chr in intersect(unique(as.character(seqnames(gr))), seqlevels(gr_snv))) {
      snv1 <- gr_snv[seqnames(gr_snv) == chr & !is.na(gr_snv$phase) & gr_snv$phase == "phase1"]
      snv2 <- gr_snv[seqnames(gr_snv) == chr & !is.na(gr_snv$phase) & gr_snv$phase == "phase2"]

      if (length(snv1) > 0)
         kpSegments(kp, chr = chr, x0 = start(snv1), x1 = start(snv1),
                    y0 = 0, y1 = 1, r0 = R0_S1, r1 = R1_S1, col = "grey20", lwd = 0.6)
      if (length(snv2) > 0)
         kpSegments(kp, chr = chr, x0 = start(snv2), x1 = start(snv2),
                    y0 = 0, y1 = 1, r0 = R0_S2, r1 = R1_S2, col = "grey20", lwd = 0.6)
   }

   # (d) Legende ancestry
   legend("bottom", legend = names(cols_ancestries), fill = cols_ancestries,
          border = NA, bty = "n", horiz = TRUE, cex = 0.7, inset = -0.02, xpd = TRUE)
}