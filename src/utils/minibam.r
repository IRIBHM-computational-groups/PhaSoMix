library(data.table)
library(tidyr)
library(parallel)

minibam_path = "rawdata/PCAWG/minibam/"
bam_dir <- "rawdata/PCAWG/bam/"
bed_dir <- "rawdata/PCAWG/bedfile_snv_mnv/"

dir.create(minibam_path, showWarnings = FALSE, recursive = TRUE)
dir.create(bed_dir,      showWarnings = FALSE, recursive = TRUE)

# Début du log
log_file <- file.path(minibam_path,"minibam_generation_log.txt")
log_con <- file(log_file, open = "w")  # ouvrir la connexion en mode append
sink(log_con, type = "message")        # rediriger les messages d'erreur
message("\n===== LOG STARTED at ", Sys.time(), " =====\n")

data_release <- fread("rawdata/PCAWG/info/pcawg-data-releases.tsv", sep = "\t", header = TRUE, data.table = FALSE)
data_release <- separate_rows(data_release, tumor_wgs_aliquot_id, tumor_wgs_bwa_alignment_bam_file_name, sep = ",")
data_release <- data_release[, c("dcc_project_code", "icgc_donor_id", "tumor_wgs_aliquot_id",
                                 "normal_wgs_bwa_alignment_bam_file_name", "tumor_wgs_bwa_alignment_bam_file_name")]

bam_files <- list.files(
  path = bam_dir,
  pattern = "\\.bam$",
  full.names = TRUE
)


mclapply(bam_files, function(bam_file_path) {
    
    message(Sys.time(), " >>> Processing ", bam_file_path)
    
    
    ID = data_release$tumor_wgs_aliquot_id[match(basename(bam_file_path),data_release$tumor_wgs_bwa_alignment_bam_file_name)]
    
    if (is.na(ID)) {
      message(Sys.time(), "  !! ID not found for ", bam_file_path)
      next
    }

    if (!file.exists(paste0(bam_file_path, ".bai"))) {
        message(Sys.time(), "  !! No index found for ", bam_file_path, ", indexing...")
        system(paste0("samtools index ", bam_file_path))
    }
    
    read_ids <- file.path(
      minibam_path,
      gsub("\\.bam$", "_read_ids.txt", basename(bam_file_path))
    )    
    
    minibam_unpaired <- file.path(
        minibam_path, 
        paste0("minibam_unpaired.PCAWG", gsub("^PCAWG\\.", "", basename(bam_file_path)))
    )
    
    minibam_paired <- file.path(
        minibam_path, 
        paste0("minibam.PCAWG.", gsub("^PCAWG\\.", "", basename(bam_file_path)))
    )
    
    if (file.exists(minibam_paired)) {
        message(Sys.time(), "  -> Minibam files already exist for ", ID)
      next
    }
    
    snv_path1 <- paste0("/srv/home/mlef0011/rawdata/PCAWG_data/TCGA/snv_mnv/",
            ID, ".consensus.20160830.somatic.snv_mnv.vcf.gz")
    snv_path2 <- paste0("/srv/home/mlef0011/rawdata/PCAWG_data/ICGC/snv_mnv/",
            ID, ".consensus.20160830.somatic.snv_mnv.vcf.gz")

    snv_path <- if (file.exists(snv_path1)) snv_path1 else if (file.exists(snv_path2)) snv_path2 else NA
   
    if (is.na(snv_path)) {
      message(Sys.time(), "  !! SNV file not found for ", ID)
      next
    }
    
    if (file.exists(snv_path)){
        
        bed_path <- file.path(
          bed_dir,
          gsub("\\.vcf\\.gz$", ".bed", basename(snv_path))
        )  
        
        # VCF vers BED
        message(Sys.time(), "  -> Converting VCF to BED")
        system(paste0("bcftools query -f '%CHROM\\t%POS\\t%POS\\n' ", 
                    snv_path," | ",
                    "awk '{print $1\"\\t\"$2-1\"\\t\"$3}' > ", bed_path))
        message(Sys.time(), "  -> VCF to BED done")

        # Génération du minibam

        message(Sys.time(), "  -> Generating minibam")
        system(paste0("samtools view -b -L ", bed_path, " ", bam_file_path, " > ", minibam_unpaired))
        system(paste0("samtools view ", minibam_unpaired, " | cut -f1 | sort | uniq > ", read_ids))

        system(paste0("samtools view -b -N ", read_ids, " ", bam_file_path, " > ", minibam_paired))
        message(Sys.time(), "  -> Minibam generated")

        # Indexation minibam
        message(Sys.time(), "  -> Indexing minibam")
        system(paste0("samtools index ", minibam_paired))
        message(Sys.time(), "  -> Minibam indexing done")
        

        # Nettoyage
        message(Sys.time(), "  -> Cleaning up temporary files")
        file.remove(minibam_unpaired)
        file.remove(read_ids)
        message(Sys.time(), "  -> Done with ", ID)
        
    } else {
                message(Sys.time(), "  !! SNV file missing for ", ID)
    }

}, mc.cores=10)
message("\n===== LOG ENDED at ", Sys.time(), " =====\n")
sink(type = "message")
close(log_con)