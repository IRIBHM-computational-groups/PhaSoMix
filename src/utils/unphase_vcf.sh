#!/bin/bash

# ============================================================
# Title:     Unphase a Phased VCF
# Author:    Maxime Lefebvre
# Created:   2025-11-26
# Purpose:   Convert a phased VCF (genotypes encoded with '|'
#            separator, e.g. 0|1) to an unphased VCF (genotypes
#            encoded with '/' separator, e.g. 0/1).
#            This step is required before re-phasing with BEAGLE5:
#            BEAGLE5 expects unphased input and will refuse to
#            process a VCF that already contains phased genotypes.
# Usage:     bash unphase_vcf.sh <input.vcf.gz>
# Inputs:    - input.vcf.gz : phased, bgzip-compressed VCF
# Outputs:   - <dir>/<prefix>_unphased.vcf.gz in the same directory
#              as the input file
# Depends:   bcftools, bgzip (htslib)
# Notes:     - The replacement is performed genome-wide with sed;
#              the VCF header is preserved unchanged.
#            - The output is bgzip-compressed for downstream
#              compatibility with BEAGLE5 and bcftools.
# ============================================================

set -euo pipefail   # exit on error, unbound variable, or pipeline failure

## ============================================================
## Argument Validation
## ============================================================

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 input.vcf.gz" >&2
    exit 1
fi

INPUT="$1"

if [ ! -f "$INPUT" ]; then
    echo "Error: file '$INPUT' not found." >&2
    exit 1
fi

## Check required tools are available
for tool in bcftools bgzip; do
    if ! command -v "$tool" &> /dev/null; then
        echo "Error: '$tool' is not installed or not in PATH." >&2
        exit 1
    fi
done

## ============================================================
## Build Output Path
## ============================================================

DIR=$(dirname  "$INPUT")
BASE=$(basename "$INPUT" .vcf.gz)
OUTPUT="${DIR}/${BASE}_unphased.vcf.gz"

echo "  -> Input  : $INPUT"
echo "  -> Output : $OUTPUT"

## ============================================================
## Unphasing Pipeline
##   bcftools view  : decompress and stream the VCF
##   sed 's/|/\//g' : replace all phased separators with unphased
##   bgzip -c        : recompress on stdout, write to OUTPUT
## ============================================================

bcftools view "$INPUT" \
    | sed 's/|/\//g' \
    | bgzip -c > "$OUTPUT"

echo "[OK] Unphased VCF written to: $OUTPUT"