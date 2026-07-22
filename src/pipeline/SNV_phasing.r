#!/usr/bin/env Rscript

# ============================================================
# Title:     SNV phasing using nearby heterozygous SNPs
# Author:    Maxime Lefebvre & Mathieu Parmentier
# Created:   01/12/2025
# Purpose:   Phase SNVs using mpileup read overlap with Het SNPs
# Requires:  R packages: data.table, GenomicRanges, stringr, parallel, dplyr
# Input:     BAM, Beagle (phased SNPs), SNV list, reference genome
# Output:    Phased SNV file per sample
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(stringr)
  library(parallel)
  library(dplyr)  
})

# Extract heterozygous SNP positions from Beagle output
get_hz_SNP <- function(chrs, beagle, sample, output_path_sample, rewrite=FALSE){
    message("   >> Extracting heterozygous SNPs")
    
      
    cols = c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", sample)
    
    output_path = file.path(output_path_sample, paste0("SNPhz_pos_", sample, ".txt"))
    if (!file.exists(output_path) | rewrite == TRUE){

        hz_snp = lapply(paste0("chr",chrs), function(chr) {
            input_file_beagle = gsub("chrCHR", chr, beagle)
            SNP_sample <- fread(input_file_beagle, skip = "#CHROM", sep = "\t", header=T, select = cols, data.table=FALSE)
            SNP_sample_hz <- subset(SNP_sample, SNP_sample[[ncol(SNP_sample)]] %in% c("0|1", "1|0"))
            pos <- SNP_sample_hz[, c(1,2,4,5)]
            return(pos)
        })
        hz_snp_df = as.data.frame(do.call(rbind,hz_snp))
        write.table(hz_snp_df, output_path, sep = "\t", col.names = FALSE, row.names = FALSE, quote=FALSE)
    } else {
        message("   -> File already exist:",output_path)
    }
    return(output_path)
}



# Run samtools mpileup on selected loci
mpileup <- function(ref_genome,bam_file, loci, output){
          
    cmd = paste("samtools mpileup",
                "-f", ref_genome,
                bam_file,
                "-Q 0",
                "-l", loci,
                "--output-QNAME -o", output, sep=" ")
    
    message("   -> Command: ", cmd)

    system(cmd, wait=T)
    return(output)
}


# Run mpileup for heterozygous SNPs
mpileup_SNP <- function(sample, bam_file, ref_genome, SNP_pos, output_path_sample, rewrite = TRUE){
    
    message("   >> Mpileup for SNP")

    output_path = file.path(output_path_sample, paste0("mplieup_SNP_",sample,".txt"))
    if (!file.exists(output_path) | rewrite == TRUE){
        mpileup(ref_genome, bam_file, SNP_pos, output_path)
        message("   -> Output written: ", output_path)
    } else {
        message("   -> File already exist:",output_path)
    }
    return(output_path)
}

# Run mpileup for SNVs
mpileup_SNV <- function(sample, bam_file, input_snv, ref_genome, output_path_sample, rewrite = TRUE){
    
    message("   >> Mpileup for SNV")

      
    output_path = file.path(output_path_sample, paste0("mplieup_SNV_",sample,".txt"))
    snv_pos_path = file.path(output_path_sample,paste0("SNV_pos_",sample,".txt"))
    
    if (!file.exists(output_path) | rewrite == TRUE){
        snv = fread(input_snv, sep = '\t', header=TRUE)
        snv = snv[snv$FILTER == ".",]
        snv_pos = snv[,c("#CHROM", "POS")]
        write.table(snv_pos,file = snv_pos_path ,sep="\t",col.names=F,row.names=F,quote=F)
        mpileup(ref_genome, bam_file, snv_pos_path, output_path)
        message("   -> Output written: ", output_path)
    } else {
        message("   -> File already exist:",output_path)
    }
    return(output_path)

}

# Filter SNPs in SNV Â± 820 bp window
FilterSNPbyWindow = function(snv, snp_tab, window = 820) {
  snv = snv[!duplicated(snv$label),]
  gr_snv = GRanges(unlist(snv$chr), IRanges(as.numeric(unlist(snv$pos)) - window, as.numeric(unlist(snv$pos)) + window)) 
  gr_snp = GRanges(snp_tab$V1, IRanges(snp_tab$V2, snp_tab$V2)) 
  overlap = gr_snp[gr_snp %over% gr_snv,]
  overlap = overlap[overlap %outside% GRanges(unlist(snv$chr), IRanges(as.numeric(unlist(snv$pos)), as.numeric(unlist(snv$pos)) ))]
  snp_filter =snp_tab[which(paste0(snp_tab$V1,":", snp_tab$V2) %in% paste0(overlap)),]
  return(snp_filter)
}

# Check if allele counts match a heterozygous profile
isHet <- function(c, pbin=0.01, probHom=.99)
{
  x= c[1] #ref
  y= c[2] #alt
  !any(pbinom(x,size=x+y,prob=probHom,log.p=F)>pbin)
}

# Parse mpileup read-base string
# Removes indels and special characters from mpileup bases.
GetParsingOfReadBase = function(ch){
  if (grepl("\\$",ch)) {
    ch = gsub("\\$", "", ch)
  }
  if (grepl("\\^",ch)) {
    ch = gsub("\\^.", "", ch)
  }
  if (grepl("\\d", ch)) {
    numbers = unlist(stringr::str_extract_all(ch, "\\d+"))
    indel = list()
    for (i in length(numbers):1) {
      n = numbers[[i]]
      indice = unlist(gregexpr(n, ch))
      if (substr(ch, indice-1, indice-1) == "+") {
        indel = append(paste0("+",substr(ch, indice-2,indice-2),substr(ch, indice + nchar(n), indice + nchar(n) - 1 + as.numeric(n))),indel)
        ch = paste0(substr(ch, 1,indice-3),"X", substr(ch, indice + nchar(n) + as.numeric(n) ,nchar(ch)))
      } else if (substr(ch, indice-1, indice-1) == "-") {
        indel = append(paste0("-",substr(ch, indice-2,indice-2),substr(ch, indice + nchar(n), indice + nchar(n) - 1 + as.numeric(n))),indel)
        ch = paste0(substr(ch, 1,indice-3),"X", substr(ch, indice + nchar(n) + as.numeric(n) ,nchar(ch)))
      }
    }
  } else {
    indel = list()
  }
  
  return(list(read_base= ch,indels = indel))
}
  
# Parse mpileup output for SNVs
# Extracts allele and read names from SNV mpileup.
GetPilueTable_SNV = function(pileup, bypass_IsHomo = FALSE) { ## work only with SNV
  
    SN = data.frame("chr"=NA, "pos" = NA, "ref" = FALSE ,"al" = NA,"al_qn" = NA, "label" = NA)
    
    nucleo = c("A","C","G","T","star") 
    allelecount = c(0,0,0,0,0) 
    Qnames_list = unlist(str_split(pileup[7], ",", n = Inf, simplify = FALSE))
    if (grepl("\\.",pileup[5])) {
      pileup[5] = gsub("\\.", pileup[3], pileup[5])
    }
    if (grepl("\\,",pileup[5])) {
      pileup[5] = gsub("\\,", pileup[3], pileup[5])
    }
    pileup[5] = toupper(pileup[5])
    reads_base = GetParsingOfReadBase(pileup[5])
    if (nchar(reads_base$read_base) != as.numeric(pileup[4])) {
      warning(paste0("\nERROR: in parsing read base (number of read base is not equal to the depth ",pileup[5]," => ",pileup[1],":",trimws(pileup[2])))
    }
    ## Do alleleCount-like
    # for single nucleotide

    al_list = unlist(str_split(reads_base$read_base, "", n = Inf, simplify = FALSE))
    allelecount[1] = length(na.omit(Qnames_list[al_list == "A"]))
    allelecount[2] = length(na.omit(Qnames_list[al_list == "C"]))
    allelecount[3] = length(na.omit(Qnames_list[al_list == "G"]))
    allelecount[4] = length(na.omit(Qnames_list[al_list == "T"]))
    allelecount[5] = length(na.omit(Qnames_list[al_list == "*"]))


    if (sum(allelecount) - sum(sort(allelecount,decreasing =TRUE)[c(1:2)]) > 3 | (sort(allelecount,decreasing =TRUE)[1] == sort(allelecount,decreasing =TRUE)[2]  & sort(allelecount,decreasing =TRUE)[2] == sort(allelecount,decreasing =TRUE)[3])) {
        
        warning(paste0("\nSKIP: more than 2 alleles detected for position => ",pileup[1],":",trimws(pileup[2])))

    } else {

        both_al = sort(allelecount,decreasing =TRUE)[c(1:2)]
        if (isHet(both_al) | bypass_IsHomo ) {

            SN = lapply(1:sum(both_al > 0), function(n) {
                chr = pileup[1]
                pos = trimws(pileup[2])
                index = which(allelecount == both_al[n])
                if (length(index) == 2 | ( n == 2 & length(index) > 2)) {
                    index = which(allelecount == both_al[n])
                    ref = nucleo[index[n]] == pileup[3]  
                    al =nucleo[index[n]]
                    al_qn = paste0(unlist(Qnames_list[al_list == nucleo[index[n]]]),collapse=";")
                } else  {
                    ref = nucleo[allelecount == both_al[n]] == pileup[3]  
                    al = nucleo[allelecount == both_al[n]]
                    al_qn = paste0(unlist(Qnames_list[al_list == nucleo[allelecount == both_al[n]]]),collapse=";")
                }
                label = paste0(pileup[1],":",trimws(pileup[2]))
                SN_a = c(chr, pos, ref, al, al_qn, label)
                return(SN_a)
            })
        SN = as.data.frame(do.call(rbind, SN))
        colnames(SN) = c("chr", "pos", "ref", "al", "al_qn", "label")
        SN[SN$al == "*"]$al = "star"

        }
    }

  
  return(SN) 
}



# Parse mpileup output for SNPs
# Extracts allele and read names from SNP mpileup.
GetPilueTable_SNP = function(pileup) {
  
    SN = data.frame("chr"=NA, "pos" = NA, "ref" = FALSE ,"al" = NA,"al_qn" = NA, "label" = NA)
    
    nucleo = c("A","C","G","T","star") 
    allelecount = c(0,0,0,0,0) 
    Qnames_list = unlist(str_split(pileup[7], ",", n = Inf, simplify = FALSE))
    if (grepl("\\.",pileup[5])) {
      pileup[5] = gsub("\\.", pileup[3], pileup[5])
    }
    if (grepl("\\,",pileup[5])) {
      pileup[5] = gsub("\\,", pileup[3], pileup[5])
    }
    pileup[5] = toupper(pileup[5])
    reads_base = GetParsingOfReadBase(pileup[5])
    if (nchar(reads_base$read_base) != as.numeric(pileup[4])) {
      warning(paste0("\nERROR: in parsing read base (number of read base is not equal to the depth ",pileup[5]," => ",pileup[1],":",trimws(pileup[2])))
    }
    ## Do alleleCount-like
    # for single nucleotide

    al_list = unlist(str_split(reads_base$read_base, "", n = Inf, simplify = FALSE))
    allelecount[1] = length(na.omit(Qnames_list[al_list == "A"]))
    allelecount[2] = length(na.omit(Qnames_list[al_list == "C"]))
    allelecount[3] = length(na.omit(Qnames_list[al_list == "G"]))
    allelecount[4] = length(na.omit(Qnames_list[al_list == "T"]))
    allelecount[5] = length(na.omit(Qnames_list[al_list == "*"]))


    if (sum(allelecount) - sum(sort(allelecount,decreasing =TRUE)[c(1:2)]) > 3 | (sort(allelecount,decreasing =TRUE)[1] == sort(allelecount,decreasing =TRUE)[2]  & sort(allelecount,decreasing =TRUE)[2] == sort(allelecount,decreasing =TRUE)[3])) {
        warning(paste0("\nSKIP: more than 2 alleles detected for position => ",pileup[1],":",trimws(pileup[2])))
        
    } else if (sum(allelecount) == 0  ){
        warning(paste0("\nSKIP: No alleles detected for position => ",pileup[1],":",trimws(pileup[2]))) 

    } else {

        both_al = sort(allelecount,decreasing =TRUE)[c(1:2)]

        SN = lapply(1:sum(both_al >0), function(n) {
            chr = pileup[1]
            pos = trimws(pileup[2])
            index = which(allelecount == both_al[n])
            if (length(index) == 2 | ( n == 2 & length(index) > 2)) {
                index = which(allelecount == both_al[n])
                ref = nucleo[index[n]] == pileup[3]  
                al =nucleo[index[n]]
                al_qn = paste0(unlist(Qnames_list[al_list == nucleo[index[n]]]),collapse=";")
            } else  {
                ref = nucleo[allelecount == both_al[n]] == pileup[3]  
                al = nucleo[allelecount == both_al[n]]
                al_qn = paste0(unlist(Qnames_list[al_list == nucleo[allelecount == both_al[n]]]),collapse=";")
            }
            label = paste0(pileup[1],":",trimws(pileup[2]))
            SN_a = c(chr, pos, ref, al, al_qn, label)
            return(SN_a)
        })
        SN = as.data.frame(do.call(rbind, SN))
        colnames(SN) = c("chr", "pos", "ref", "al", "al_qn", "label")
        SN[SN$al == "*"]$al = "star"

    }
  
  return(SN) 
}


# Infer SNVâ€“SNP phasing from shared reads
# Determines phase by overlapping read names.
GetPhasing = function(snv,snp) {
  

    SNV = snv[which(snv$ref == FALSE),]
    #SNP = snp[which(snp$ref == FALSE),]
    SNP <- unique(snp[,c("chr","pos")])
    window = 820  
    gr_snv = GRanges(unlist(SNV$chr), IRanges(as.numeric(unlist(SNV$pos)) - window, as.numeric(unlist(SNV$pos)) + window)) 
    gr_snp = GRanges(unlist(SNP$chr), IRanges(as.numeric(unlist(SNP$pos)), as.numeric(unlist(SNP$pos))))

    phase = mclapply(1:length(SNV$chr), function(i) {
        flag_overlap = gr_snp %over% gr_snv[c(i),]
        cand_snv = snv[which(snv$chr == SNV$chr[i] & snv$pos == SNV$pos[i]),]
        Qname_snv_al   = unlist(str_split(cand_snv[which(cand_snv$ref == FALSE),5],";", n = Inf, simplify = FALSE))
        Qname_snv_ref  = unlist(str_split(cand_snv[which(cand_snv$ref == TRUE ),5],";", n = Inf, simplify = FALSE))

        if (sum(flag_overlap) > 0){ 

            overlap = gr_snp[flag_overlap,]
            phase_SNP = lapply(1:length(overlap), function(j) {
    
                phase_a = NULL

                index = which(snp$chr == as.character(seqnames(overlap[j])) & snp$pos == as.character(ranges(overlap[j])))
                cand_snp = snp[index,]
                if (nrow(cand_snp) > 0){
                    Qname_snp_al  = unlist(str_split(cand_snp[which(cand_snp$ref == FALSE),5],";", n = Inf, simplify = FALSE)) 
                    Qname_snp_ref = unlist(str_split(cand_snp[which(cand_snp$ref == TRUE ),5],";", n = Inf, simplify = FALSE)) 

                    SNPref_SNVref = sum(Qname_snv_ref %in% Qname_snp_ref) 
                    SNPalt_SNVref = sum(Qname_snv_ref %in% Qname_snp_al)
                    SNPref_SNValt = sum(Qname_snv_al %in% Qname_snp_ref)
                    SNPalt_SNValt = sum(Qname_snv_al %in% Qname_snp_al)

                    if (sum(SNPref_SNValt, SNPalt_SNValt) > 0) {

                        phase_a = data.frame("chr_snv"=cand_snv[which(cand_snv$ref == FALSE),]$chr, 
                                             "pos_snv" = cand_snv[which(cand_snv$ref == FALSE),]$pos,
                                             "reads_snv" = length(Qname_snv_al), 
                                             "chr_snp"= unique(cand_snp$chr), 
                                             "pos_snp" = unique(cand_snp$pos),  
                                             "SNPref_SNVref" = sum(Qname_snv_ref %in% Qname_snp_ref), 
                                             "SNPalt_SNVref" = sum(Qname_snv_ref %in% Qname_snp_al),
                                             "SNPref_SNValt" = sum(Qname_snv_al %in% Qname_snp_ref), 
                                             "SNPalt_SNValt" = sum(Qname_snv_al %in% Qname_snp_al))
                    }
                }

                return(phase_a)   
            })
            
            phase_SNP = Filter(Negate(is.null), phase_SNP)
            if (length(phase_SNP) > 0) {
              phase_SNP = do.call(rbind, phase_SNP)
            } else {
              phase_SNP = data.frame()
            }
            return(phase_SNP)
        } else {
            return(NULL)
        }
    },mc.cores = 10)
    
    phase = Filter(Negate(is.null), phase)
    if (length(phase) > 0) {
      phase = do.call(rbind, phase)
    } else {
      phase = data.frame()
    }
 
    phase$SNPref_SNVref <- as.numeric(phase$SNPref_SNVref)
    phase$SNPalt_SNVref <- as.numeric(phase$SNPalt_SNVref)
    phase$SNPref_SNValt <- as.numeric(phase$SNPref_SNValt)
    phase$SNPalt_SNValt <- as.numeric(phase$SNPalt_SNValt)

    phase$phasing = apply(phase, 1, get_phasing)
    #phase <- phase[!is.na(phase$phasing), ]
    phase$SNV_with_SNP <- phase$SNPalt_SNValt > phase$SNPref_SNValt
    
    phase_output <- phase %>%
        distinct() %>%
        group_by(chr_snv, pos_snv, reads_snv) %>%
        summarise(across(everything(), ~ paste(., collapse = ","))) %>%
        ungroup()
    phase_output = as.data.frame(phase_output)
  return(phase_output)
}

# Assign final phase to each SNV
get_phasing <- function(row) {
    SNV_state <- c(as.numeric(row["SNPalt_SNValt"]), as.numeric(row["SNPref_SNValt"]))

    if (sum(SNV_state) > 0) {
        p_val <- suppressWarnings(wilcox.test(c(rep(1, as.numeric(SNV_state[1])), rep(-1, as.numeric(SNV_state[2]))))$p.value)
        if ((SNV_state[1] == 0 & SNV_state[2] > 0) | 
            (SNV_state[2] == 0 & SNV_state[1] > 0) | 
            (SNV_state[1] <= 3 & p_val < 0.05) | 
            (SNV_state[2] <= 3 & p_val < 0.05)) {
            phasing <- paste(SNV_state, collapse = "|")
        } else {
            phasing <- NA
        }
    } else {
        phasing <- NA
    }
    
    return(phasing)
}



# Run full SNV phasing pipeline for one sample
# Runs mpileup, parsing, filtering and phasing.
phase_SNV <- function(chrs, sample, input_bam, input_vcf, output_path_sample, input_SNV, ref_genome_path, rewrite=FALSE){
    
    message(">>> SNV phasing for sample:", sample)

    phased_SNV_path = file.path(output_path_sample,paste0('SNV_phasing_mpileup_',sample,".txt"))
    if(file.exists(phased_SNV_path) & !rewrite){
        message("   -> File already exist:",phased_SNV_path)
        return(phased_SNV_path)
    }
    
    # Check BAM and reference
    if (!file.exists(input_bam)) stop("Missing BAM file: ", input_bam)
    if (!file.exists(ref_genome_path)) stop("Missing reference genome: ", ref_genome_path)

    # Check Beagle and SNV files for chr1â€“23
    missing <- unlist(lapply(chrs, function(chr) {
        beagle_chr <- sub("chrCHR", paste0("chr", chr), input_vcf)
        snv_chr    <- sub("chrCHR", paste0("chr", chr), input_SNV)
        c(if (!file.exists(beagle_chr)) beagle_chr, if (!file.exists(snv_chr)) snv_chr)
    }))

    if (length(missing) > 0) {
        stop("Missing files:\n", paste(missing, collapse = "\n"))
    }

    cat("âœ“ All files found. Proceeding...\n")
                    
    if (!dir.exists(output_path_sample)) dir.create(output_path_sample, recursive = TRUE)

    # mpileup
    SNP_pos = get_hz_SNP(chrs, input_vcf, sample, output_path_sample)
    input_mpileup_snp = mpileup_SNP(sample, input_bam, ref_genome_path, SNP_pos, output_path_sample, rewrite = rewrite)
    input_mpileup_snv = mpileup_SNV(sample, input_bam, input_SNV, ref_genome_path, output_path_sample,rewrite = rewrite)


    # mpileup parsing
    message("  >> SNV Mpileup output parsing" )
    pileups_snv = fread(input_mpileup_snv, sep = "\t")
    pileups_snv$V1 = as.character(pileups_snv$V1)
    SNV <- apply(na.omit(pileups_snv), 1, function(row) GetPilueTable_SNV(row))
    SNV = as.data.frame(do.call(rbind, SNV))
    SNV = na.omit(SNV)
    
    message("  >> SNP Mpileup output parsing" )
    pileups_snp = fread(input_mpileup_snp, sep = "\t")
    pileups_snp$V1 = as.character(pileups_snp$V1)
    pileups_snp_filter= FilterSNPbyWindow(SNV,pileups_snp)
    SNP<- apply(na.omit(pileups_snp_filter), 1, function(row) GetPilueTable_SNP(row))
    SNP = as.data.frame(do.call(rbind,SNP))
    SNP = na.omit(SNP)

    # SNV phasing
    message("  >> SNVs phasing with SNPs" )
    phase = GetPhasing (SNV,SNP)
    phase <- phase[order(as.numeric(phase$chr_snv)), ]

    # write output
    write.table(phase, phased_SNV_path, sep = "\t", col.names = TRUE, row.names = FALSE, quote=FALSE)
    message("   -> Output written: ", phased_SNV_path)

    message(">>> SNV phasing for sample ", sample, " completed.")
    return(phased_SNV_path)
}