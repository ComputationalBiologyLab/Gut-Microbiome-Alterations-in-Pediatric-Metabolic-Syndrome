# Gut microbiome alterations in pediatric metabolic syndrome: harmonized multi-cohort reproducibility analysis


Three public pediatric fecal 16S rRNA datasets (PRJNA1016449, USA; PRJNA838778, China; PRJNA819279, Mexico) are reprocessed with one pipeline, tested for differential abundance with ANCOM-BC2, pooled by random-effects meta-analysis and evaluated with a pre-specified reproducibility classification. Raw reads are available from the NCBI Sequence Read Archive under the BioProject accessions above. Derived outputs are in the Supplementary Tables of the article.

## Repository layout

```
environment/            conda environment (QIIME 2 amplicon 2024.10, DADA2, cutadapt, R 4.3.3 and packages)
pipeline/<BioProject>/  single-cohort pipeline
    config/             analysis_config.yaml (primers, region, platform, group definitions), env.sh
    scripts/            numbered steps 00-08 and R scripts
meta_analysis/          pooling, reproducibility classification, phylum-level analysis, figures
sensitivity_analyses/   robustness, power, concordance benchmark, alternative estimators,
                        attenuation outside the discovery cohort, rarefaction sensitivity
```

## Requirements

- conda environment: `conda env create -n qiime2-amplicon-2024.10 -f environment/qiime2-amplicon-2024.10.yml`
- Python 3.10 or later with the packages in `requirements.txt` (meta-analysis, figures and sensitivity analyses)
- R 4.3.3 with phyloseq 1.46.0, vegan 2.6.8 and ANCOMBC 2.4.0 (included in the conda environment)
- SILVA 138 99% naive Bayes classifier for QIIME 2 (downloaded by `05_taxonomy.sh` when absent)

## Running the analysis

Each cohort has its own working directory named `<BioProject>_analysis` with the folders `00_metadata/` (sample metadata and read manifest) and `01_raw_data/` (FASTQ files downloaded from the SRA).

1. Single-cohort pipeline. From the cohort directory, with the scripts and config folders of `pipeline/<BioProject>/` copied into it:
   ```bash
   bash scripts/00_run_pipeline.sh
   ```
   Steps: `01_verify_raw_data.sh` (checksums and manifest), `02_raw_qc.sh`, `03_preprocess_cutadapt.sh` (primer removal), `04_dada2.sh` (denoising), `05_taxonomy.sh` (SILVA 138), `06_filter_and_phylogeny.sh` (ASV filters, tree), `07_diversity.sh` (rarefaction, diversity), `08_downstream_r.sh` (genus aggregation, ANCOM-BC2, effect-size tables).
2. Meta-analysis. Set `MS_ROOT` to the directory that contains the three `<BioProject>_analysis` folders:
   ```bash
   export MS_ROOT=/path/to/project
   python meta_analysis/run_meta_analysis.py            # pooling and reproducibility classification
   Rscript meta_analysis/run_phylum_ancombc2_all_studies.R
   python meta_analysis/run_phylum_meta_analysis.py     # phylum-level sensitivity analysis
   python meta_analysis/make_prisma_chart.py
   python meta_analysis/make_manuscript_figures.py
   ```
3. Sensitivity analyses (read the pipeline and meta-analysis outputs through `MS_ROOT`):
   ```bash
   cd sensitivity_analyses
   python robustness.py        # pairwise agreement, Hartung-Knapp, fixed effect, leave-one-cohort-out
   python power.py             # minimum detectable effect, power, equivalence intervals
   python concordance_null.py  # assumptions behind the 25% chance benchmark
   python reestimate.py        # alternative tau2 estimators, Q-profile intervals, Bayesian model
   python attenuation.py       # reduction of effects outside the discovery cohort
   Rscript rarefaction_sensitivity.R
   ```

## Key parameters

- Primer removal: cutadapt through QIIME 2 (error rate 0.15, indels and wildcards allowed, minimum length 50). Primers: 338F/806R (PRJNA1016449), 515F/806R (PRJNA838778, single-end analysis), 341F/805R (PRJNA819279).
- DADA2: independent pooling, consensus chimera removal, maximum expected errors 2, truncation quality 2; paired cohorts truncated at 260/200 bp; PRJNA838778 not truncated by length.
- Taxonomy: SILVA 138 99% naive Bayes classifier, confidence 0.7; mitochondria and chloroplast sequences removed.
- ASV filter: present in at least ceil(5% of samples) samples and at least 10 reads in total. Genus prevalence of at least 10% within a cohort.
- Rarefaction (diversity only): 20th percentile of per-sample read counts (at least 5,000 reads), replaced by the 10th percentile (at least 1,000) if fewer than 80% of samples would be kept.
- Differential abundance: ANCOM-BC2 on unrarefied genus counts, model `~ mets_status` (US sensitivity model `~ mets_status + age_years + sex + BMI`), BH adjustment, `struc_zero = TRUE`, `neg_lb = TRUE`, `s0_perc = 0.05`, `alpha = 0.05`, `prv_cut = 0`, `lib_cut = 0`, `pseudo_sens = FALSE`. Significance within a cohort: FDR < 0.1.
- Meta-analysis: genera present in all cohorts with finite estimates; DerSimonian-Laird random effects, BH across genera. Classes: strongly concordant (same direction in all cohorts, FDR < 0.20, I2 < 60%), directionally concordant, heterogeneous, no clear association (all absolute effects below 0.25 and p >= 0.05).
- Random seeds are set in the scripts (for example 1016449 for PRJNA1016449).

## Notes

- The pipeline scripts were written for macOS and use the conda environment prefix defined in `pipeline/<BioProject>/config/env.sh` (override with `MINICONDA_ROOT`).
- Group labels follow the primary studies: `Control` in the code denotes the non-MetS comparator group of each cohort.

- `pipeline/PRJNA1016449/scripts/downstream_analysis.R` also contains an optional MaAsLin2 comparison that runs only when the MaAsLin2 package is installed; it was not used for any result in the article (ANCOM-BC2 is the only differential-abundance method reported).

## License

MIT (see `LICENSE`).
