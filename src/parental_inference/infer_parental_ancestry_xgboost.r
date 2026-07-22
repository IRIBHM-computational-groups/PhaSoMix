setwd("/srv/home/mlef0011/Phasomix/")
source("src/utils/XGBOOST.r")
library(parallel)

info_ma = read.table("/srv/home/mlef0011/Phasomix/output/local_ancestry/prediction/high_coverage/IS_corrected/input_path_F1_like_mono_x_admixed_gnomix_PCAWG_POP_no_AMR.tsv", header=T, sep = "\t")

input_vcf_sex = "rawdata/PCAWG/rawdata/chrCHR_PCAWG_snp_subset_1kgp_PCAWG_admixed_no_AMR.vcf.gz"
input_model_xgboost ="output/parental_origin/training/POP_no_AMR/"


samples <- unique(info_ma$sample)
parental_origin_inference_list <- lapply(samples, function(sample){
    hapA = sort(unlist(strsplit(gsub("[{}]", "", info_ma[info_ma$sample == sample, 'hapA']), ",")))
    hapB = sort(unlist(strsplit(gsub("[{}]", "", info_ma[info_ma$sample == sample, 'hapB']), ",")))
    sexe <- info_ma[info_ma$sample == sample, "sexe"]


    ancestry_mtDNA <- get_xgboost_prediction_prob(
        sample, input_vcf_sex, "chrMT", input_model_xgboost
    )

    

    ancestry_y <- get_xgboost_prediction_prob(
        sample, input_vcf_sex, "chrY", input_model_xgboost
    )

    

    ancestry_x <- get_xgboost_prediction_prob(
        sample, input_vcf_sex, "chrX", input_model_xgboost
    )
    
    return(list(sample = sample,
                sexe = sexe,
                hapA = hapA,
                hapB = hapB,
                ancestry_mtDNA = ancestry_mtDNA,
                ancestry_y = ancestry_y,
                ancestry_x = ancestry_x))
})