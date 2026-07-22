# PhaSoMix

**Parent-of-origin phasing of somatic mutations in admixed cancer genomes**

PhaSoMix assigns somatic mutations to their parental haplotype of origin in admixed cancer patients, without requiring parent or parent-surrogate sequencing. It exploits the genetic divergence between maternal and paternal haplotypes in recently admixed individuals as a proxy for parent-of-origin, and is applied here to the Pan-Cancer Analysis of Whole Genomes (PCAWG) cohort to test whether the two parental genomes accumulate somatic mutations at equal rates.

> Lefèbvre M, Cléris A, Parmentier M, Van Loo P, Detours V, Tarabichi M. *Parent-of-origin phasing of somatic mutations shows equal mutation burden between parental genomes in human cancers.* (manuscript in preparation)

---

## Overview

In admixed individuals, the two parental haplotypes carry distinct local ancestries across part of the genome. PhaSoMix uses these ancestry-discriminant regions as a stand-in for parental identity, and combines local ancestry inference, SNP-in-read phasing, and sex-chromosome/mtDNA ancestry classification to phase somatic mutations to the maternal or paternal genome. Accuracy is benchmarked and propagated at every step to bound any true difference in mutation burden between the two parental haplotypes.

The pipeline is organized around five tasks:

1. **Global ancestry inference** (ADMIXTURE) — identify ancestry-pure individuals to enrich the local ancestry reference panel.
2. **SNP phasing and switch-error correction** (Beagle 5.2 + allelic-imbalance-based correction).
3. **Local ancestry inference** (Gnomix) and identification of informative (admixed) patients via bipartite graph resolution of parental haplotypes.
4. **Parental origin inference** (XGBoost on mtDNA/chrX/chrY) to resolve which haplotype is maternal and which is paternal.
5. **Somatic mutation phasing and quantification**, corrected for allele-specific copy number, with per-sample and per-region outputs.

---

## Repository structure

This reflects the current content of the repository:

```
PhaSoMix/
├── input/
│   ├── config_0.1cM.yaml                              # Gnomix configuration (window size, smoothing, calibration, …)
│   └── input_path_F1_like_mono_x_admixed_gnomix_PCAWG_POP_no_AMR.tsv
│                                                        # Sample sheet: sample, ancestry, project, hapA/hapB,
│                                                        # sexe, paths to BAM/allelecount/SNV/CNA/CCF files
│
└── src/
    ├── local_ancestry/
    │   ├── training_gnomix.r                # Train Gnomix local ancestry models (super-population pairs)
    │   ├── training_gnomix_sex_chr.r        # Train Gnomix models restricted to sex-chromosome analyses
    │   ├── training_gnomix_sim.r            # Train Gnomix models for simulation benchmarks
    │   ├── prediction_gnomix_PCAWG.r        # Predict local ancestry on admixed PCAWG patients
    │   ├── prediction_gnomix_sim.r          # Predict local ancestry on simulated individuals
    │   ├── accuracy_gnomix_sim.r            # Local ancestry accuracy vs. switch-error rate (simulations)
    │   └── get_informative_patients.r       # Identify informative patients; parental haplotype resolution
    │                                          via co-occurrence graph + bipartite colouring (igraph);
    │                                          karyoplots of local ancestry per patient
    │
    ├── parental_inference/
    │   ├── training_XGboost_parental_origin.r     # Train XGBoost models per ancestry pair on mtDNA/chrX/chrY
    │   ├── infer_parental_ancestry_xgboost.r       # Resolve maternal/paternal haplotype for each patient
    │   └── accuracy_XGboost_parental_origin.r      # Classification accuracy on held-out 1kGP samples
    │
    ├── pipeline/
    │   ├── SNV_quantification.r             # Main CLI driver — runs the 4 steps below per sample
    │   ├── SNV_phasing.r                    # Phase somatic SNVs to nearby heterozygous SNPs (mpileup overlap)
    │   ├── SNV_ancestry.r                   # Assign each phased SNV to a local-ancestry haplotype
    │   └── SNV_correction.r                 # Determine gained-copy status and correct SNV multiplicity for CN
    │
    ├── utils/
    │   ├── utility_functions.r              # Shared I/O, validation, and ancestry helper functions
    │   ├── writevcf.r                       # VCF writer (header + bgzip compression)
    │   ├── ADMIXTURE.r                      # PLINK/ADMIXTURE wrapper functions
    │   ├── score_SNPs.r                     # Ancestry-informativeness scoring of SNPs from 1kGP allele frequencies
    │   ├── extract_discriminant_SNPs.r      # Select the most ancestry-discriminant SNPs
    │   ├── XGBOOST.r                        # XGBoost training/prediction helpers for chrX/chrY/mtDNA
    │   ├── genotyping_sex_chr.R             # Genotype mtDNA/chrX/chrY from 1kGP VCFs for XGBoost training
    │   ├── allele_count_sex_chr.r           # Run alleleCounter on PCAWG BAMs at sex-chromosome/mtDNA loci
    │   ├── minibam.r                        # Extract minibams around SNV/MNV positions for SNV phasing
    │   ├── compute_SER_IS.r                 # Switch-error rate computation in imbalanced segments (PCF + LLR)
    │   ├── correct_SE_IS.r                  # Apply switch-error correction across imbalanced segments
    │   ├── get_BAF.r                        # B-allele frequency computation from allele counts
    │   ├── get_stats.r                      # Per-sample SNV counts by ancestry / CNV zone
    │   ├── local_quantification.r           # Region-level quantification (genes/exons/introns/imprinted/CpG)
    │   ├── get_circos_plot.r                # Circos plot of SNV ancestry + CNV tracks per sample
    │   ├── get_karyoplot.r                  # Karyoplot of Gnomix local ancestry calls per sample
    │   ├── simulations.r                    # Simulate mixed-ancestry genomes with controlled switch-error rates
    │   ├── phase_beagle.sh                  # Beagle 5.2 SNP phasing against 1kGP reference panels, per chromosome
    │   └── unphase_vcf.sh                   # Strip phasing (`|` → `/`) from a VCF before re-phasing
    │
    └── accuracy/                            # Jupyter notebooks reproducing the accuracy/power benchmarks
        ├── accuracy_local_ancestry.ipynb        #   Local ancestry inference accuracy (Gnomix, per ancestry pair/SER)
        ├── accuracy_SNV_phasing.ipynb           #   SNV-to-SNP phasing accuracy
        ├── accuracy_sex_chr.ipynb               #   mtDNA/chrX/chrY parental-origin classification accuracy
        ├── accuracy_correction.ipynb            #   Copy-number correction accuracy (simulated + PCAWG)
        ├── SE_rate_imbalance_zones.ipynb        #   Switch-error rate before/after correction (low/high coverage)
        ├── comparison_pred_SNV_corrected_vs_bo_corrected.ipynb  # SE-correction impact on ancestry/SNV assignment
        └── statictical_power.ipynb              #   Statistical power / null-distribution analysis
```

> `rawdata/`, `output/`, and `software/` are used throughout the scripts as relative paths (1kGP/PCAWG reference data, pipeline outputs, and third-party binaries such as PLINK, ADMIXTURE, and Beagle) but are not tracked in this repository, since they contain controlled-access data and large external binaries. See [Requirements](#requirements) and [Data availability](#data-availability) below to reconstitute them locally.

---

## Requirements

### External tools
- [Beagle 5.2](https://faculty.washington.edu/browning/beagle/beagle.html) — SNP phasing
- [Gnomix](https://github.com/AI-sandbox/gnomix) — local ancestry inference
- [ADMIXTURE](https://dalexander.github.io/admixture/) — global ancestry inference
- [PLINK 1.9](https://www.cog-genomics.org/plink/)
- [bcftools](https://samtools.github.io/bcftools/) / [samtools](http://www.htslib.org/) / `bgzip` (htslib)
- [alleleCounter](https://github.com/cancerit/alleleCount)

### R (≥ 4.x)
```r
install.packages(c("data.table", "dplyr", "tidyr", "stringr", "parallel",
                    "igraph", "optparse", "RColorBrewer", "circlize"))

# Bioconductor
BiocManager::install(c("GenomicRanges", "Biostrings", "karyoploteR"))

install.packages("xgboost")
```

### Python / Jupyter
The notebooks in `src/accuracy/` require a Jupyter kernel with `R` (via `IRkernel`) or Python, depending on the notebook; see each notebook's first cell for its specific dependencies.

---

## Pipeline workflow

**1. Global ancestry inference**
```r
source("src/utils/score_SNPs.r")
source("src/utils/extract_discriminant_SNPs.r")
source("src/utils/ADMIXTURE.r")
```
Scores 1kGP SNPs by ancestry informativeness, extracts the most discriminant subset, and runs ADMIXTURE to flag ancestry-pure PCAWG individuals (dominant ancestry > 0.95), which are added to the Gnomix training set.

**2. SNP phasing and switch-error correction**
```bash
bash src/utils/phase_beagle.sh
```
```r
source("src/utils/correct_SE_IS.r")   # calls compute_SER_IS.r / get_BAF.r internally
```
Phases 1kGP + PCAWG genotypes with Beagle 5.2, then detects and corrects residual switch errors within allelic-imbalance segments using BAF-based piecewise-constant fitting and a likelihood-ratio test.

**3. Local ancestry inference and informative patient identification**
```r
source("src/local_ancestry/training_gnomix.r")
source("src/local_ancestry/prediction_gnomix_PCAWG.r")
source("src/local_ancestry/get_informative_patients.r")
```
Trains Gnomix per super-population pair, predicts local ancestry along each phased PCAWG haplotype, and identifies informative (admixed) patients by resolving parental haplotypes through a co-occurrence graph and bipartite colouring.

**4. Parental origin inference**
```r
source("src/parental_inference/training_XGboost_parental_origin.r")
source("src/parental_inference/infer_parental_ancestry_xgboost.r")
```
Trains XGBoost classifiers on mtDNA/chrX/chrY genotypes per ancestry pair, then resolves which haplotype (hapA/hapB) is maternal and which is paternal for each informative patient.

**5. Somatic mutation phasing and quantification**
```bash
Rscript src/pipeline/SNV_quantification.r \
    -i input/input_path_F1_like_mono_x_admixed_gnomix_PCAWG_POP_no_AMR.tsv \
    -o output/ \
    -v "rawdata/vcf/chrCHR_phased.vcf.gz" \
    -s "rawdata/vcf/chrCHR_sex_phased.vcf.gz" \
    -r rawdata/reference/hg19.fa \
    -g output/local_ancestry/prediction/ \
    -x output/parental_origin/training/ \
    -t 10
```
The main driver. For each sample it: (1) phases somatic SNVs to nearby heterozygous SNPs, (2) assigns each phased SNV to a local-ancestry haplotype, (3) corrects SNV counts for allele-specific copy number, (4) generates a circos plot, and (5) computes global and region-level (gene/exon/intron/imprinted/CpG) mutation statistics per parental haplotype.

| Flag | Description |
|---|---|
| `-i, --path_info_ma` | Sample sheet TSV (see `input/`) |
| `-o, --output_path` | Output directory |
| `-v, --vcf_ma` | Phased autosomal VCF, `CHR` placeholder for chromosome |
| `-s, --vcf_ma_sex` | Phased sex-chromosome/mtDNA VCF |
| `-r, --ref_genome_path` | Reference genome FASTA (hg19) |
| `-g, --gnomix_pred_path` | Gnomix prediction directory |
| `-x, --xgboost_model_path` | XGBoost parental-origin model directory |
| `-a, --accuracy` | Accuracy summary file path |
| `-t, --threads` | Number of threads (default: 10) |

**Accuracy benchmarks.** `src/local_ancestry/accuracy_gnomix_sim.r`, `src/parental_inference/accuracy_XGboost_parental_origin.r`, and the notebooks in `src/accuracy/` reproduce the simulation- and PCAWG-based benchmarks (local ancestry accuracy vs. switch-error rate, SNV phasing accuracy, sex-chromosome/mtDNA classification, copy-number correction, and statistical power) used to bound the detectable mutation-rate asymmetry between parental genomes.

---

## Input data

- **`input/config_0.1cM.yaml`** — Gnomix configuration (0.1 cM window size, smoothing, context ratio, calibration flag).
- **`input/input_path_F1_like_mono_x_admixed_gnomix_PCAWG_POP_no_AMR.tsv`** — sample sheet consumed by `SNV_quantification.r`, with one row per sample: `sample`, `ancestry`, `project`, `hapA`/`hapB` (ancestry sets per haplotype), `sexe`, and paths to the BAM, alleleCounter output, somatic SNV VCF, copy-number segments, and subclonal/CCF file.

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

Maxime Lefèbvre — maxime.lefebvre@ulb.be  
Maxime Tarabichi — maxime.tarabichi@ulb.be

Issues and questions are welcome via the [GitHub issue tracker](https://github.com/IRIBHM-computational-groups/PhaSoMix/issues).
