# PhaSoMix

**Parent-of-origin phasing of somatic mutations in admixed cancer genomes**

PhaSoMix is a computational pipeline that assigns somatic mutations to their parental haplotype of origin in admixed cancer patients, without requiring parent or parent-surrogate sequencing. It exploits the genetic divergence between maternal and paternal haplotypes in recently admixed individuals as a proxy for parent-of-origin, and applies it to the Pan-Cancer Analysis of Whole Genomes (PCAWG) cohort to test whether the two parental genomes accumulate somatic mutations at equal rates.

This repository contains the full pipeline described in:

> Lefèbvre M, Cléris A, Parmentier M, Van Loo P, Detours V, Tarabichi M. *Parent-of-origin phasing of somatic mutations shows equal mutation burden between parental genomes in human cancers.* (manuscript in preparation)

---

## Overview

In admixed individuals, the two parental haplotypes carry distinct local ancestries across part of the genome. PhaSoMix uses these ancestry-discriminant regions as a stand-in for parental identity, and combines local ancestry inference, SNP-in-read phasing, and sex-chromosome/mtDNA ancestry classification to phase somatic mutations to the maternal or paternal genome. Errors are quantified and propagated at every step to provide an upper bound on any true difference in mutation burden between the two parental haplotypes.

The pipeline performs four main tasks:

1. **Local ancestry inference** across the cohort to identify informative (admixed) patients, i.e. those in whom at least one ancestry uniquely tags a single parental haplotype.
2. **Empirical phasing of somatic mutations** to nearby heterozygous SNPs using raw sequencing reads.
3. **Parental (maternal/paternal) origin assignment** from mitochondrial DNA and chromosome X/Y ancestry.
4. **Aggregation of phased mutations** by parental origin into genome-wide mutation burdens, corrected for allele-specific copy number.

Simulation-based accuracy benchmarks (local ancestry inference, SNV phasing, sex-chromosome/mtDNA classification, copy-number correction) are used to propagate uncertainty through the pipeline and bound the detectable mutation-rate asymmetry between parental genomes.

---

## Repository structure

```
PhaSoMix/
├── src/
│   ├── admixed_inference.r                    # Identify admixed individuals from ADMIXTURE output
│   ├── admixed_inference_accuracy.r            # Repeated CV accuracy of admixture-based classification
│   ├── training_gnomix_paired.r                # Train Gnomix local ancestry models (pop/subpop pairs)
│   ├── training_gnomix_paired_sim.r            # Train Gnomix models for simulation benchmarks
│   ├── prediction_gnomix_paired_mixed_ancestry.r  # Local ancestry prediction on admixed PCAWG patients
│   ├── prediction_gnomix_paired_simulations.r  # Local ancestry prediction on simulated individuals
│   ├── accuracy_gnomix_paired.r                # Gnomix accuracy across ancestry pairs / switch-error rates
│   ├── accuracy_sex_chr_paired.r               # XGBoost accuracy for chrX/chrY/mtDNA ancestry classification
│   ├── training_XGboost_parental_origin.r      # Train XGBoost models for mtDNA/chrX/chrY ancestry
│   ├── correct_SE_IS.r                         # Switch-error detection and correction in imbalanced segments
│   ├── simulations_mixed_ancestry.r            # Simulate admixed genomes with controlled switch-error rates
│   ├── SNV_quantification.r                    # Main CLI driver: phasing → ancestry → parental origin → quantification
│   └── utils/
│       ├── utility_functions.r                 # Shared I/O, validation, and ancestry helper functions
│       ├── writevcf.r                          # VCF writer (header + compression)
│       ├── ADMIXTURE.r                         # PLINK/ADMIXTURE wrapper functions
│       ├── extract_discriminant_SNPs.r         # Ancestry-informative SNP selection
│       ├── XGBOOST.r                           # XGBoost training/prediction for sex chromosomes and mtDNA
│       ├── infer_parental_origin.r             # Resolve maternal/paternal haplotype from mtDNA/chrX/chrY
│       ├── SNV_phasing.r                       # SNV-to-SNP read-based phasing (mpileup overlap)
│       ├── SNV_ancestry.r                      # Assign phased SNVs to a local-ancestry haplotype
│       ├── get_BAF.r                           # B-allele frequency computation from allele counts
│       ├── compute_SER_IS.r                    # Switch-error rate computation in imbalanced segments (PCF + LLR)
│       ├── SNV_correction.r                    # Copy-number correction of SNV multiplicity
│       ├── get_stats.r                         # Per-sample SNV counts by ancestry/CNV zone
│       ├── local_quantification.r              # Region-level (gene/exon/intron/imprinted/CpG) quantification
│       └── quantification_plotting.r           # Circos and summary plots per sample
├── utils/
│   └── unphase_vcf.sh                          # Strip phasing (| → /) from a VCF before re-phasing
├── phase_PCAWG_mixed_ancestry.sh               # Beagle 5.2 phasing of admixed PCAWG samples, per chromosome
├── software/                                   # Third-party binaries (PLINK, ADMIXTURE, Beagle .jar) — not versioned
├── rawdata/                                    # Reference panels, genetic maps, annotation files (1kGP, PCAWG)
├── input/                                      # Sample sheets, Gnomix config files
└── output/                                     # Pipeline results (VCFs, ancestry calls, quantification tables, figures)
```

> `software/`, `rawdata/`, `input/`, and `output/` are not tracked in version control (data and binaries); see [Requirements](#requirements) and [Data availability](#data-availability) for how to populate them.

---

## Requirements

### External tools
- [Beagle 5.2](https://faculty.washington.edu/browning/beagle/beagle.html) (SNP phasing)
- [Gnomix](https://github.com/AI-sandbox/gnomix) (local ancestry inference)
- [ADMIXTURE](https://dalexander.github.io/admixture/) (global ancestry inference)
- [PLINK 1.9](https://www.cog-genomics.org/plink/)
- [bcftools](https://samtools.github.io/bcftools/) / [samtools](http://www.htslib.org/)
- [alleleCounter](https://github.com/cancerit/alleleCount)
- [pyega3](https://github.com/EGA-archive/ega-download-client) (EGA data download)

### R (≥ 4.x)
```r
install.packages(c("data.table", "dplyr", "tidyr", "stringr", "parallel",
                    "ggplot2", "patchwork", "ggnewscale", "circlize",
                    "reshape2", "RColorBrewer", "optparse", "igraph"))

# Bioconductor
BiocManager::install(c("GenomicRanges", "Biostrings", "karyoploteR", "copynumber"))

install.packages("xgboost")
install.packages("lme4")  # downstream mutation-asymmetry testing
```

---

## Pipeline workflow

The pipeline is organized into five stages, run roughly in this order:

**1. Global ancestry & informative patient identification**
```bash
Rscript src/admixed_inference.r
```
Runs ADMIXTURE on discriminant SNPs and identifies admixed PCAWG patients whose parental haplotypes carry distinguishable ancestries (strict F1 or partially informative).

**2. Phasing and switch-error correction**
```bash
bash phase_PCAWG_mixed_ancestry.sh
Rscript src/correct_SE_IS.r
```
Phases genotypes with Beagle 5.2 against the 1kGP high-coverage reference panel, then detects and corrects residual switch errors within allelic-imbalance segments using BAF-based piecewise-constant fitting and a likelihood-ratio test.

**3. Local ancestry inference**
```bash
Rscript src/training_gnomix_paired.r
Rscript src/prediction_gnomix_paired_mixed_ancestry.r
```
Trains Gnomix models per super-population pair on 1kGP + ancestry-pure PCAWG haplotypes, then predicts local ancestry along each phased haplotype of the admixed patients.

**4. Parental origin resolution**
```bash
Rscript src/training_XGboost_parental_origin.r
```
Trains XGBoost classifiers on mtDNA, chrX, and chrY genotypes to resolve maternal versus paternal ancestry for each informative patient (called internally via `infer_parental_origin()` in the next step).

**5. Somatic mutation phasing and quantification**
```bash
Rscript src/SNV_quantification.r \
    -i path_info_mixed_ancestry.tsv \
    -o output/ \
    -v "rawdata/vcf/chrCHR_phased.vcf.gz" \
    -r rawdata/reference/hg19.fa
```
The main driver script: phases somatic SNVs to nearby heterozygous SNPs (read-based), assigns each phased SNV to a local-ancestry haplotype, resolves parental origin, corrects mutation counts for copy number, and produces per-sample and per-region (gene/exon/intron/imprinted/CpG) quantification tables and plots.

Accuracy and simulation scripts (`accuracy_gnomix_paired.r`, `accuracy_sex_chr_paired.r`, `admixed_inference_accuracy.r`, `simulations_mixed_ancestry.r`, `prediction_gnomix_paired_simulations.r`) reproduce the benchmarking experiments used to quantify and propagate error at each pipeline step.

---

## Data availability

- **1000 Genomes Project (1kGP):** low-coverage Phase 3 (GRCh37) and high-coverage (GRCh38) panels, publicly available at [ftp.1000genomes.ebi.ac.uk](https://ftp.1000genomes.ebi.ac.uk/).
- **PCAWG:** whole-genome BAM files are controlled-access, available through the [European Genome-phenome Archive](https://ega-archive.org/) (EGAS00001001692) for ICGC samples and [Bionimbus/ICGC](https://icgc.bionimbus.org/files) for TCGA-derived samples. Somatic variant calls are available via the [ICGC ARGO platform](https://platform.icgc-argo.org/) and Bionimbus/ICGC.

Access to controlled data requires separate data access approval through ICGC/TCGA.

---

## Citation

If you use PhaSoMix, please cite:

> Lefèbvre M, Cléris A, Parmentier M, Van Loo P, Detours V, Tarabichi M. Parent-of-origin phasing of somatic mutations shows equal mutation burden between parental genomes in human cancers. (manuscript in preparation)

---

## Contact

Maxime Lefèbvre — IRIBHM, Université Libre de Bruxelles (ULB)
Maxime Tarabichi — maxime.tarabichi@ulb.be

Issues and questions are welcome via the [GitHub issue tracker](https://github.com/IRIBHM-computational-groups/PhaSoMix/issues).
