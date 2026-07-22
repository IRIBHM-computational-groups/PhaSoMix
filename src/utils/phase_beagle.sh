#!/bin/bash
MEM="100g"

for chr in {1..22}; do

    ref_file_lc="rawdata/1kGp/rawdata/chr${chr}.1kGp_low_coverage_subset_PCAWG_snp_unrelated_IDs_1kGp.vcf.gz"
    ref_file_hc="rawdata/1kGp/rawdata/chr${chr}.1kGp_high_coverage_Illumina_subset_1kgp_PCAWG_snp_unrelated_IDs_1kGp.vcf.gz"
    #ref_file_1kGp_PCAWG="rawdata/ref_files/ref_files_1kGp_high_coverage_PCAWG_no_MA_corrected_ISs_chr${chr}.vcf.gz"

    map_file="/srv/home/mlef0011/Phasomix/rawdata/1kGp/mapfiles/plink.chr${chr}.GRCh37.map"

    gt="rawdata/PCAWG/rawdata/chr${chr}_PCAWG_beagle5_subset_1kgp_PCAWG_snp.vcf.gz"

    echo ">>> Traitement chr${chr}..."

    # ---------- UNPHASE ----------
    bash src/utils/unphase_vcf.sh "$gt"
    gt_unphased="rawdata/PCAWG/rawdata/chr${chr}_PCAWG_beagle5_subset_1kgp_PCAWG_snp.vcf.gz"

    # ---------- OUTPUT NAMES ----------
    base=${gt%.vcf.gz}
    output_lc="${base}_phased_1kGp_unrelared_IDs_low_coverage_reference_panel"
    output_hc="${base}_phased_1kGp_unrelared_IDs_high_coverage_reference_panel"
    #output_1kGp_PCAWG="${base}_1kGp_high_coverage_corrected_PCAWG_reference_panel"

    # ---------- BEAGLE - low coverage ----------
    java -Xmx${MEM} -Xms4g -XX:+UseParallelOldGC -jar software/beagle.22Jul22.46e.jar \
        gt="$gt_unphased" \
        ref="$ref_file_lc" \
        out="$output_lc" \
        map="$map_file" \
        nthreads=10 \
        window=40 \
        overlap=4 \
        impute=true

    # ---------- BEAGLE - high coverage ----------
    java -Xmx${MEM} -Xms4g -XX:+UseParallelOldGC -jar software/beagle.22Jul22.46e.jar \
        gt="$gt_unphased" \
        ref="$ref_file_hc" \
        out="$output_hc" \
        map="$map_file" \
        nthreads=10 \
        window=40 \
        overlap=4 \
        impute=true

    # ---------- BEAGLE - 1kGp high coverage + PCAWG corrected ----------
    #java -Xmx${MEM} -Xms4g -XX:+UseParallelOldGC -jar software/beagle.22Jul22.46e.jar \
     #   gt="$gt_unphased" \
     #   ref="$ref_file_1kGp_PCAWG" \
     #   out="$output_1kGp_PCAWG" \
     #   map="$map_file" \
     #   nthreads=4 \
     #   window=40 \
     #   overlap=4 \
     #  impute=true

    if [ $? -ne 0 ]; then
        echo "ERREUR sur chr${chr}" >&2
        continue
    fi

    echo ">>> chr${chr} terminé."
done