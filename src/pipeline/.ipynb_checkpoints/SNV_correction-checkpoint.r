# ============================================================
# Title:     SNV Quantification with Copy Number Correction
# Author:    Maxime Lefebvre
# Created:   01/12/2025
# Purpose:   Determine SNV gain status and compute corrected allelic multiplicity
#            based on ancestry, Beagle phasing, BAF, and copy number data
# Requires:  R packages: data.table
# Output:    Corrected SNV quantification files per sample
# ============================================================

source("src/utils/get_BAF.r")

# Determine if a SNV is on the gained copy
is_gain <- function(chr, sample, input_beagle, input_allelecount,IZ_SNVs, output_path){
    
    cols = c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", sample)

    input_file_beagle = gsub("chrCHR", paste0("chr",chr), input_beagle)
    input_file_allelecount = gsub("chrCHR", paste0("chr",chr), input_allelecount)
    beagle = fread(input_file_beagle, skip = "#CHROM", sep = "\t", header=T, select = cols, data.table=FALSE)
    BAF_raw = get_baf(sample, chr, input_file_allelecount, beagle, output_path = output_path)
    BAF = merge(beagle,BAF_raw, by.x = c('#CHROM','POS'), by.y = c("#CHR","POS"))
    IZ_SNVs_chr =  IZ_SNVs[IZ_SNVs$chr_snv == as.character(chr),] 


    if (nrow(IZ_SNVs_chr) > 0){

        IZ_SNVs_chr$gain_copy = sapply(1:nrow(IZ_SNVs_chr), function(i) {

            chr = as.vector(unlist(strsplit(IZ_SNVs_chr[i,"chr_snp"],split = ",")))
            pos_snp = as.vector(unlist(strsplit(IZ_SNVs_chr[i,"pos_snp"],split = ",")))
            SNV_with_SNP = as.vector(unlist(strsplit(IZ_SNVs_chr[i,"SNV_with_SNP"],split = ",")))

            gain_copy = sapply(1:length(chr), function(idx){

                chr_split = chr[idx]
                pos_snp_split = pos_snp[idx]
                SNV_with_SNP_split = SNV_with_SNP[idx]

                phase = BAF[BAF$"#CHROM" == chr_split & BAF$POS == pos_snp_split,10]
                BAF_phase1 = BAF[BAF$"#CHROM" == chr_split & BAF$POS == pos_snp_split,"BAF_phase1"]
                BAF_phase2 = BAF[BAF$"#CHROM" == chr_split & BAF$POS == pos_snp_split,"BAF_phase2"]

                gain_copy = FALSE

                if (SNV_with_SNP_split == TRUE){
                    if (phase == "0|1"){
                        if (BAF_phase2 > BAF_phase1){
                            gain_copy = TRUE
                        }
                    } else if (phase == "1|0"){
                        if (BAF_phase1 > BAF_phase2){
                            gain_copy = TRUE
                        }
                    }
                } else if (SNV_with_SNP_split == FALSE){
                    if (phase == "1|0"){
                        if (BAF_phase2 > BAF_phase1){
                            gain_copy = TRUE
                        }
                    } else if (phase == "0|1"){
                        if (BAF_phase1 > BAF_phase2){
                            gain_copy = TRUE
                        }
                    }
                }
                return(gain_copy)
            })
            return(ifelse(all(gain_copy == gain_copy[1]), gain_copy[1], NA))
        })
        return(IZ_SNVs_chr)
    } else {
        return(NULL)
    }

}
                                                                                                                    
# Correct SNV quantification based on copy number
get_correction  <- function(IZ_SNVs){
    IZ_SNVs$correction = ifelse(
        IZ_SNVs$gain_copy,
        IZ_SNVs$mult/IZ_SNVs$major_cn, 
        ifelse(
            IZ_SNVs$minor_cn == 0,
            0,
            IZ_SNVs$mult/IZ_SNVs$minor_cn
        )
    )
    return(IZ_SNVs)
}

#Quantify SNV allelic multiplicity genome-wide                                                                               
get_quantification <- function(sample, ancestry_SNV_path, input_SNV_ccf, input_beagle, input_allelecount, output_path, rewrite=FALSE){
    
    message(">>> SNVs quantification for sample:", sample)
   
    quantification_SNV_path = file.path(output_path,paste0('SNV_quantification_',sample,".txt"))

    if(file.exists(quantification_SNV_path) & !rewrite){
        message("   -> File already exist:",quantification_SNV_path)
        return(quantification_SNV_path)
    }  
    
    # Check files
    if (!file.exists(ancestry_SNV_path)) stop("Missing ancestry SNV file: ", ancestry_SNV_path)
    if (!file.exists(input_SNV_ccf)) stop("Missing SNV info file: ", input_SNV_ccf)
    
    SNVs = fread(ancestry_SNV_path, header=T, sep= "\t", data.table=F)
    SNVs <- subset(SNVs, !is.na(ancestry_SNV))
    SNVs$chr_snv = as.character(SNVs$chr_snv)

    # Check Beagle and SNV files for chr1–22
    missing <- unlist(lapply(1:22, function(chr) {
        beagle_chr <- sub("chrCHR", paste0("chr", chr), input_beagle)
        allelecount_chr    <- sub("chrCHR", paste0("chr", chr), input_allelecount)
        c(if (!file.exists(beagle_chr)) beagle_chr,
        if (!file.exists(allelecount_chr)) allelecount_chr)
    }))

    common_chrs <- unlist(lapply(1:22, function(chr) {
        beagle_chr <- sub("chrCHR", paste0("chr", chr), input_beagle)
        allelecount_chr    <- sub("chrCHR", paste0("chr", chr), input_allelecount)
        if (file.exists(beagle_chr) && file.exists(allelecount_chr)) chr
    }))

    if (length(missing) > 0) {
        cat("Missing files:\n", paste(missing, collapse = "\n"), "\n")
    }

    if (length(common_chrs) == 0) stop("No common chromosomes found")

    cnp <- fread(input_SNV_ccf,sep = '\t', header=TRUE, data.table=F)

    SNVs_cnp = merge(cnp[,c("chromosome","position","major_cn","minor_cn","mult")],SNVs , by.x = c("chromosome","position"), by.y = c("chr_snv","pos_snv"))
        
    colnames(SNVs_cnp)[1:2] = colnames(SNVs)[1:2]

    no_IZ_SNVs = subset(SNVs_cnp, !is.na(SNVs_cnp$minor_cn) & !is.na(SNVs_cnp$major_cn) & SNVs_cnp$minor_cn == SNVs_cnp$major_cn)
    no_IZ_SNVs$type = "no_IZ"
    no_IZ_SNVs$gain_copy = NA
    no_IZ_SNVs$correction <- no_IZ_SNVs$mult / no_IZ_SNVs$minor_cn


    IZ_SNVs = subset(SNVs_cnp, !is.na(SNVs_cnp$minor_cn) & !is.na(SNVs_cnp$major_cn) & SNVs_cnp$minor_cn != SNVs_cnp$major_cn)
    IZ_SNVs$type = ifelse(IZ_SNVs$minor_cn == 0, "IZ_LOH","IZ_no_LOH")

    # Determination if the mutation is on the gain copy in imbalance zones
    list_IZ_SNVs_gain_copy = lapply(intersect(common_chrs,unique(SNVs$chr_snv)), function(chr) is_gain(chr, sample, input_beagle, input_allelecount ,IZ_SNVs, output_path))
    list_IZ_SNVs_gain_copy <- list_IZ_SNVs_gain_copy[!sapply(list_IZ_SNVs_gain_copy, is.null)]

    if (length(list_IZ_SNVs_gain_copy) > 0) {

        IZ_SNVs <- do.call(rbind, list_IZ_SNVs_gain_copy)
        # Correction of mutation quantification regarding CN 
        message("  >> Correction for CN")
        IZ_SNVs = get_correction(IZ_SNVs)

    } else { 

        IZ_SNVs = data.frame()
    }

    # Join data (IZ and no_IZ)
    SNVs_quant = rbind(no_IZ_SNVs, IZ_SNVs)
    SNVs_quant$chr_snv = as.numeric(SNVs_quant$chr_snv)
    SNVs_quant <- SNVs_quant[order(SNVs_quant$chr_snv, SNVs_quant$pos_snv ), ]
    write.table(SNVs_quant, quantification_SNV_path, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
    message(">>> SNV quantification for sample ",sample, "is done.")
    return(quantification_SNV_path)
}