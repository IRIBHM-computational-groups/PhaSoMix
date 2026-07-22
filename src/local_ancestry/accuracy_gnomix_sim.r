library(dplyr)
library(tidyr)
library(data.table)
library(parallel)

source("src/utils/utility_functions.r")

get_SE_genotype <- function(n_snps, se_pos, sample){

   ancestries_vec <- strsplit(sub("_.*", "", sample), "-")[[1]]
   ancestry1 <- ancestries_vec[1]
   ancestry2 <- ancestries_vec[2]

   true_gt <- rep(paste0(ancestry1, "|", ancestry2), n_snps)
   starts <- sort(se_pos[se_pos$INDV == sample, "IDX"])

   if (length(starts) == 0) return(true_gt)

   ends <- c(starts[-1] - 1, n_snps)
   flip_intervals <- seq_along(starts) %% 2 == 1

   mask <- logical(n_snps)
   for (k in which(flip_intervals)) {
      mask[starts[k]:ends[k]] <- TRUE
   }

   SE_gt <- true_gt
   SE_gt[mask] <- paste0(ancestry2, "|", ancestry1)

   return(SE_gt)
}

get_predicted_genotype <- function(sample, pred_label, pop_vector) {

   hap_cols <- c(
      paste0(sample, ".0"),
      paste0(sample, ".1")
   )

   n_snps <- pred_label$`n snps`

   haplotype <- paste(
      pop_vector[pred_label[[hap_cols[1]]] + 1],
      pop_vector[pred_label[[hap_cols[2]]] + 1],
      sep = "|"
   )

   return(rep(haplotype, n_snps))
}


OUTPUT_DIR    <- "output/local_ancestry/accuracy/prediction/1kGp_only/"
PRED_PATH     <- "output/local_ancestry/accuracy/prediction/1kGp_only/"
POP_INFO_PATH <- "rawdata/1kGp/info/1kGP.3202_samples.pop_info.txt"
SIM_DIR       <- "input/simulations"

SE_RATES <- c(0,
  1/1000,  # 0.001  (0.1%)
  1/500,   # 0.002  (0.2%)
  1/200,   # 0.005  (0.5%)
  1/100    # 0.01   (1%)
)

SEG_LENGTH <- 0.1
CHRS <- 22:22

pop_info     <- read.table(POP_INFO_PATH,    header = TRUE, sep = "\t")

ancestries_pop <- unique(na.omit(pop_info$`Super.Population`))
ancestries_pop_no_AMR <- setdiff(unique(na.omit(pop_info$`Super.Population`)), "AMR")

list_subpop <- lapply(ancestries_pop_no_AMR, function(p) {
   unique(na.omit(pop_info$Population[pop_info$`Super.Population` == p]))
})
names(list_subpop) = ancestries_pop_no_AMR
list_ancestries <- c(list(POP_no_AMR = ancestries_pop_no_AMR),list_subpop)

mclapply(SE_RATES, function(SE){

   tryCatch({

      lapply(names(list_ancestries), function(label){

         ancestries <- list_ancestries[[label]]

         accuracy_chr_list <- lapply(CHRS, function(chr){

            tryCatch({

               message(">>> Running SE ", SE, " - Ancestry [", label, "] - CHROMOSOME [", chr, "]")

               input_pred_prob <- file.path(
                  PRED_PATH,
                  label,
                  paste0("SE_", round(SE, 5)),
                  "prediction",
                  paste0(SEG_LENGTH, "cM"),
                  paste0("chr", chr),
                  "query_results.fb"
               )

               input_pred_label <- file.path(
                  PRED_PATH,
                  label,
                  paste0("SE_", round(SE, 5)),
                  "prediction",
                  paste0(SEG_LENGTH, "cM"),
                  paste0("chr", chr),
                  "query_results.msp"
               )

               pred_prob  <- fread(input_pred_prob,  data.table = FALSE, header = TRUE)
               pred_label <- fread(input_pred_label, data.table = FALSE, header = TRUE)

               line1 <- readLines(input_pred_label, n = 1)
               line1 <- sub("^#Subpopulation order/codes:\\s*", "", line1)
               pops <- strsplit(line1, "\t")[[1]]
               pop_vector <- sub("=.*", "", pops)

               samples <- colnames(pred_label)[7:ncol(pred_label)]
               samples <- unique(sub("\\.[01]$", "", samples))
                
               n_snps <- sum(pred_label$`n snps`)

               if (SE != 0) {
                   
                   se_pos_path <- file.path(SIM_DIR,paste0("SE_positions_simulations_mixed_ancestry_",round(SE,5),"_chr",chr,".tsv"))
               se_pos <- read.table(se_pos_path, header = TRUE, sep = "\t")

                   pair_split <- data.table::tstrsplit(sub("_.*", "", se_pos$INDV), "-", fixed = TRUE)
                   keep <- pair_split[[1]] %in% pop_vector & pair_split[[2]] %in% pop_vector
                   se_pos <- se_pos[keep, ]                   
                   
                  geno_matrix_true <- do.call(
                     cbind,
                     lapply(samples, function(s) get_SE_genotype(n_snps, se_pos, s))
                  )
               } else {
                  geno_matrix_true <- do.call(
                     cbind,
                     lapply(samples, function(s){
                        ancestries_vec <- strsplit(sub("_.*", "", s), "-")[[1]]
                        rep(paste0(ancestries_vec[1], "|", ancestries_vec[2]), n_snps)
                     })
                  )
               }
               colnames(geno_matrix_true) <- samples

               geno_matrix_predicted <- do.call(
                  cbind,
                  lapply(samples, function(s) get_predicted_genotype(s, pred_label, pop_vector))
               )
               colnames(geno_matrix_predicted) <- samples

               if (ncol(geno_matrix_predicted) != ncol(geno_matrix_true) ||
                   nrow(geno_matrix_predicted) != nrow(geno_matrix_true) ||
                   !all(colnames(geno_matrix_predicted) == colnames(geno_matrix_true))) {

                  message("Dimension mismatch: SE=", SE,
                          " ancestry_lvl=", label,
                          " chr=", chr)
                  return(NULL)
               }

               accuracy_samples <- colMeans(
                  geno_matrix_predicted == geno_matrix_true,
                  na.rm = TRUE
               )

               # Renvoyer une data.frame nommée — robuste au filtrage par chr
               data.frame(
                  SE             = SE,
                  chr            = chr,
                  ancestry_level = label,
                  as.list(accuracy_samples),
                  check.names    = FALSE,
                  stringsAsFactors = FALSE
               )

            }, error = function(e) {

               message("ERROR: SE=", SE,
                       " ancestry_lvl=", label,
                       " chr=", chr,
                       " | ", conditionMessage(e))

               return(NULL)
            })

         })

         accuracy_chr_list <- Filter(Negate(is.null), accuracy_chr_list)
         if (length(accuracy_chr_list) == 0) return(invisible(NULL))

         # rbindlist gère le cas où les chr n'ont pas exactement les mêmes samples
         accuracy_chr <- rbindlist(accuracy_chr_list, fill = TRUE)

         output_file <- file.path(
            OUTPUT_DIR,
            paste0("accuracy_gnomix_simulations_mixed_ancestry_SE_",
                   round(SE, 5), "_", label, ".tsv")
         )

         fwrite(accuracy_chr, output_file, sep = "\t", quote = FALSE)

         invisible(NULL)
      })

   }, error = function(e) {

      message("FATAL ERROR in SE=", SE,
              " | ", conditionMessage(e))

      return(NULL)
   })

}, mc.cores = length(SE_RATES))