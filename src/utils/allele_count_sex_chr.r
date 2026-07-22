library(data.table)
library(tidyr)
library(parallel)

bam_dir            <- "rawdata/PCAWG/bam/"
allele_count_path  <- "rawdata/PCAWG/allelecount/"
alleleCounter_exec <- "/srv/home/mathparm/alleleCount/bin/alleleCounter"

loci_file_X  <- "rawdata/1kGp/loci/loci_1kGp_hg19_snp_non_PAR_chrX.txt"
loci_file_Y  <- "rawdata/1kGp/loci/loci_1kGp_hg19_snp_non_PAR_chrY.txt"
loci_file_MT <- "rawdata/1kGp/loci/loci_1kGp_hg19_snp_chrMT.txt"

dir.create(allele_count_path, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(allele_count_path, "allelecount_generation_log.txt")
log_con  <- file(log_file, open = "w")
sink(log_con, type = "message")
message("\n===== LOG STARTED at ", Sys.time(), " =====\n")

data_release <- fread("rawdata/PCAWG/info/pcawg-data-releases.tsv", sep = "\t", header = TRUE, data.table = FALSE)
data_release <- separate_rows(data_release, tumor_wgs_aliquot_id, tumor_wgs_bwa_alignment_bam_file_name, sep = ",")
data_release <- data_release[, c("dcc_project_code", "icgc_donor_id", "tumor_wgs_aliquot_id",
                                 "normal_wgs_bwa_alignment_bam_file_name", "tumor_wgs_bwa_alignment_bam_file_name")]

bam_files <- list.files(path = bam_dir, pattern = "\\.bam$", full.names = TRUE)

mclapply(bam_files, function(bam_file_path) {

   message(Sys.time(), " >>> Processing ", bam_file_path)

   ID <- data_release$tumor_wgs_aliquot_id[match(basename(bam_file_path),
                                                  data_release$tumor_wgs_bwa_alignment_bam_file_name)]

   if (is.na(ID)) {
      message(Sys.time(), "  !! ID not found for ", bam_file_path)
      return(NULL)
   }

   if (!file.exists(paste0(bam_file_path, ".bai"))) {
      message(Sys.time(), "  !! No index found for ", bam_file_path, ", indexing...")
      system(paste0("samtools index ", bam_file_path))
   }

   output_x  <- file.path(allele_count_path, paste0(ID, "_alleCount_chrX.txt"))
   output_y  <- file.path(allele_count_path, paste0(ID, "_alleCount_chrY.txt"))
   output_mt <- file.path(allele_count_path, paste0(ID, "_alleCount_chrMT.txt"))

   if (file.exists(output_x) && file.exists(output_y) && file.exists(output_mt)) {
      message(Sys.time(), "  -> AlleleCount files already exist for ", ID)
      return(NULL)
   }

   allele_x  <- paste(alleleCounter_exec, "-b", bam_file_path, "-o", output_x,  "-l", loci_file_X,  "-m 20 -q 35 --dense-snps")
   allele_y  <- paste(alleleCounter_exec, "-b", bam_file_path, "-o", output_y,  "-l", loci_file_Y,  "-m 20 -q 35 --dense-snps")
   allele_mt <- paste(alleleCounter_exec, "-b", bam_file_path, "-o", output_mt, "-l", loci_file_MT, "-m 20 -q 35 --dense-snps")

   message(Sys.time(), "  -> AlleleCount chrX for ", ID)
   system(allele_x)
   message(Sys.time(), "  -> AlleleCount chrY for ", ID)
   system(allele_y)
   message(Sys.time(), "  -> AlleleCount mtDNA for ", ID)
   system(allele_mt)

}, mc.cores = min(detectCores() - 1, length(bam_files)))

message("\n===== LOG ENDED at ", Sys.time(), " =====\n")
sink(type = "message")
close(log_con)