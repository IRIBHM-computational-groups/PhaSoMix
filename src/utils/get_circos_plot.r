suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(RColorBrewer)
    library(circlize)
})

ancestry_colors <- c(
   AMR = "#F0A202",
   EUR = "#378ADD",
   AFR = "#D85A30",
   EAS = "#1D9E75",
   SAS = "#7F77DD"
)

get_circos_plot <- function(sample, ancestry_SNV_path, input_snv, input_CNA, pred_path,
                            output_dir, ref_genome = "hg19", CNV_MAX = 5,
                            rewrite = FALSE, png_dpi = 600){

   output_pdf <- file.path(output_dir, paste0("circos_", sample, ".pdf"))
   output_png <- file.path(output_dir, paste0("circos_", sample, ".png"))

   if (file.exists(output_pdf) && file.exists(output_png) && !rewrite) {
      message("   -> Plots already exist: ", output_pdf, " / ", output_png)
      return(invisible(c(pdf = output_pdf, png = output_png)))
   }

   # Tous les SNV (pour le denominateur et la categorie NA)
   snv <- fread(input_snv, sep = "\t", header = TRUE)
   snv = snv[snv$FILTER == ".",]
   all_SNV <- snv[, c("#CHROM", "POS")]
   colnames(all_SNV) <- c("V1", "V2")
   all_SNV$V1 <- paste0("chr", all_SNV$V1)
   all_SNV$V3 <- all_SNV$V2 + 1

   # SNV annotes par haplotype / ancestry
   snv_tmp <- fread(ancestry_SNV_path, header = TRUE, sep = "\t", data.table = FALSE)

   snv_table_hap1 <- snv_tmp[!is.na(snv_tmp$haplotype_SNV) & snv_tmp$haplotype_SNV == "hapA", ]
   snv_table_hap1 <- data.frame(
      "V1"       = paste0("chr", snv_table_hap1$chr_snv),
      "V2"       = snv_table_hap1$pos_snv,
      "V3"       = snv_table_hap1$pos_snv + 1,
      "ancestry" = snv_table_hap1$ancestry_SNV,
      stringsAsFactors = FALSE
   )

   snv_table_hap2 <- snv_tmp[!is.na(snv_tmp$haplotype_SNV) & snv_tmp$haplotype_SNV == "hapB", ]
   snv_table_hap2 <- data.frame(
      "V1"       = paste0("chr", snv_table_hap2$chr_snv),
      "V2"       = snv_table_hap2$pos_snv,
      "V3"       = snv_table_hap2$pos_snv + 1,
      "ancestry" = snv_table_hap2$ancestry_SNV,
      stringsAsFactors = FALSE
   )

   snv_assigned <- rbind(snv_table_hap1, snv_table_hap2)
   snv_table_NA <- anti_join(all_SNV, snv_assigned, by = c("V1", "V2"))
   snv_table_NA$ancestry <- NA_character_

   # ------------------------------------------------------------
   # Decoupage par ancestry + NA
   # ------------------------------------------------------------
   ancestries_present <- sort(unique(snv_assigned$ancestry))

   SNV_table_all <- lapply(ancestries_present, function(a) {
      snv_assigned[snv_assigned$ancestry == a, c("V1", "V2", "V3")]
   })
   names(SNV_table_all) <- ancestries_present
   SNV_table_all[["NA"]] <- snv_table_NA[, c("V1", "V2", "V3")]

   clean_names <- names(SNV_table_all)
   n_total <- nrow(all_SNV)
   counts  <- sapply(SNV_table_all, nrow)

   keep          <- counts > 0
   SNV_table_all <- SNV_table_all[keep]
   clean_names   <- clean_names[keep]
   counts        <- counts[keep]

   names(SNV_table_all) <- paste0(
      clean_names, ": ", counts, " SNVs (",
      round(counts / n_total * 100, 2), "%)"
   )

   # ============================================================
   # Copy Number
   # ============================================================
   cnv <- fread(input_CNA, sep = "\t", header = TRUE, data.table = FALSE)
   cnv <- subset(cnv, !is.na(minor_cn) & !is.na(major_cn))
   cnv$chr     <- paste0("chr", cnv$chromosome)
   cnv$segmean <- cnv$major_cn + cnv$minor_cn

   cnv_cleaned <- data.frame(
      "chr"   = cnv$chr,
      "start" = cnv$start,
      "end"   = cnv$end,
      "val1"  = cnv$segmean
   )
   cnv_cleaned <- cnv_cleaned[which(cnv_cleaned$val1 <= CNV_MAX), ]
   cnv_cleaned <- cnv_cleaned[order(cnv_cleaned$val1, decreasing = FALSE), ]

   if (length(which(cnv_cleaned$val1 == 2)) == 0) {
      df_cnv2_fake <- data.frame("chr" = "chr1", "start" = 1, "end" = 2, "val1" = 2)
      cnv_cleaned  <- rbind(df_cnv2_fake, cnv_cleaned)
   }
   cnv_cleaned$legend <- paste("CN=", cnv_cleaned$val1, sep = "")

   CNV_table  <- list()
   index_list <- 1
   for (i in unique(cnv_cleaned$val1)) {
      CNV_table[[index_list]] <- cnv_cleaned[which(cnv_cleaned$val1 == i), ]
      index_list <- index_list + 1
   }

   # ============================================================
   # Ancestry locale (GNomix .msp) - phases 1 et 2
   # ============================================================
   all_chr_data <- lapply(1:22, function(chr) {
      pred_file <- file.path(pred_path, paste0("chr", chr), "query_results.msp")
      if (file.exists(pred_file)) {
         line1 <- readLines(pred_file, n = 1)
         line1 <- sub("^#Subpopulation order/codes:\\s*", "", line1)
         pops  <- strsplit(line1, "\t")[[1]]
         pop_vector <- sub("=.*", "", pops)
         cols <- c("#chm", "spos", "epos", "sgpos", "egpos", "n snps",
                   paste0(sample, c(".0", ".1")))
         pred_gnomix <- fread(pred_file, select = cols, header = TRUE, sep = "\t",
                              data.table = FALSE)
         pred_gnomix[, 7:ncol(pred_gnomix)] <- lapply(
            pred_gnomix[, 7:ncol(pred_gnomix)],
            function(col) pop_vector[col + 1]
         )
         return(pred_gnomix)
      }
   })
   all_chr_data <- Filter(Negate(is.null), all_chr_data)
   common_cols  <- Reduce(intersect, lapply(all_chr_data, colnames))
   pred_all_chr <- as.data.frame(do.call(rbind, lapply(all_chr_data,
                                                       function(df) df[, common_cols])))

                                                       
                                                       
   # palette ancestry
   pred_file <- file.path(pred_path, "chr22/query_results.msp")
   line1 <- readLines(pred_file, n = 1)
   line1 <- sub("^#Subpopulation order/codes:\\s*", "", line1)
   pops  <- strsplit(line1, "\t")[[1]]
   pops_present <- sub("=.*", "", pops)

                                                       
   cols_ancestries = ancestry_colors[pops_present]
                                                                                                      
   df <- pred_all_chr[, c(colnames(pred_all_chr)[1:6],
                          paste0(sample, c(".0", ".1")))]
   colnames(df)[(ncol(df) - 1):ncol(df)] <- c("pred_phase1", "pred_phase2")
   df$`#chm` <- paste0("chr", df$`#chm`)

   phase1_bed <- data.frame(chr   = df$`#chm`,
                            start = df$spos,
                            end   = df$epos,
                            anc   = df$pred_phase1,
                            stringsAsFactors = FALSE)
   phase2_bed <- data.frame(chr   = df$`#chm`,
                            start = df$spos,
                            end   = df$epos,
                            anc   = df$pred_phase2,
                            stringsAsFactors = FALSE)

   # ============================================================
   # Couleurs SNV
   # ============================================================
   rainfall_cols <- ifelse(
      clean_names == "NA",
      "grey",
      cols_ancestries[clean_names]
   )
   names(rainfall_cols) <- names(SNV_table_all)

   missing_anc <- setdiff(clean_names[clean_names != "NA"], names(cols_ancestries))
   if (length(missing_anc)) {
      warning("Ancestry sans couleur dans cols_ancestries : ",
              paste(missing_anc, collapse = ", "))
   }

   # ============================================================
   # FONCTION DE TRACE (appelee une fois par device)
   # ============================================================
   draw_circos <- function() {

      circos.clear()                       # reset entre deux devices
      circos.initializeWithIdeogram(species = ref_genome)
      circos.par("track.height" = 0.12)

      is_na_cat  <- grepl("^NA:", names(SNV_table_all))
      draw_order <- order(!is_na_cat)

      circos.genomicRainfall(
         SNV_table_all[draw_order],
         cex = 0.4,
         col = adjustcolor(rainfall_cols[draw_order], alpha.f = 0.55),
         pch = 16
      )

      circos.genomicDensity(
         SNV_table_all,
         col = adjustcolor(rainfall_cols, alpha.f = 0.5)
      )

      circos.par("track.margin" = c(0, 0))

      circos.genomicTrack(phase1_bed, ylim = c(0, 1), bg.border = NA, track.height = 0.07,
         panel.fun = function(region, value, ...) {
            circos.genomicRect(region, value, ytop = 1, ybottom = 0,
                               col = cols_ancestries[value$anc], border = NA, ...)
         })

      circos.genomicTrack(phase2_bed, ylim = c(0, 1), bg.border = NA, track.height = 0.07,
         panel.fun = function(region, value, ...) {
            circos.genomicRect(region, value, ytop = 1, ybottom = 0,
                               col = cols_ancestries[value$anc], border = NA, ...)
         })

      circos.par("track.margin" = c(0, 0.01))

      map <- setNames(c("orange", "black", "purple", "skyblue", "darkgreen"),
                      c("CN=1", "CN=2", "CN=3", "CN=4", "CN=5"))
      map.legend   <- map[sort(unique(unlist(lapply(CNV_table,
                                                    function(x) x %>% select(legend)))))]
      palette_cols <- as.character(map.legend)

      if (length(CNV_table) == 1) {
         circos.genomicTrack(CNV_table,
                             ylim = c(0, 1),
                             panel.fun = function(region, value, ...) {
                                i <- getI(...)
                                circos.genomicLines(region, value, col = palette_cols[i],
                                                    type = "segment", lwd = 2, ...)
                             })
      } else {
         circos.genomicTrack(CNV_table,
                             panel.fun = function(region, value, ...) {
                                i <- getI(...)
                                circos.genomicLines(region, value, col = palette_cols[i],
                                                    type = "segment", lwd = 2, ...)
                             })
      }

      legend("topleft", names(SNV_table_all), col = rainfall_cols,
             bty = "n", title = expression(bold("SNV ancestry")), cex = 0.8, pch = 19)

      legend("bottomleft", names(cols_ancestries), col = cols_ancestries,
             bty = "n", title = expression(bold("Local ancestry")), cex = 0.8, pch = 15)

      legend("topright", names(map.legend), lty = 1, lwd = 3, col = map.legend,
             bty = "n", title = expression(bold("Total Copy Number")), cex = 0.8)

      title(main = paste0("PCAWG sample:\n ", sample), cex.main = 0.8, line = -1)
   }

   # ============================================================
   # Export PDF (vectoriel - pour archive / re-edition)
   # ============================================================
   pdf(output_pdf, width = 9, height = 9)
   draw_circos()
   dev.off()

   # ============================================================
   # Export PNG haute resolution (raster - pour Inkscape)
   # 9 in x 600 dpi = 5400 px ; type = "cairo" -> anti-aliasing
   # ============================================================
   png(output_png, width = 9, height = 9, units = "in",
       res = png_dpi)
   draw_circos()
   dev.off()

   message("   -> Wrote: ", output_pdf)
   message("   -> Wrote: ", output_png, " (", png_dpi, " dpi)")

   invisible(c(pdf = output_pdf, png = output_png))
}