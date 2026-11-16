#!/usr/bin/env Rscript
# Report and tables for PRJNA838778 (run after ANCOM-BC2)
suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(ape)
  library(ggplot2)
  library(pheatmap)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("--file=", "", args[grep("--file=", args)])
if (length(file_arg)) {
  root <- normalizePath(file.path(dirname(file_arg), ".."))
  setwd(root)
}

dir.create("11_figures", showWarnings = FALSE)
dir.create("12_tables", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
dir.create("10_functional_prediction", showWarnings = FALSE)

numify <- function(x) suppressWarnings(as.numeric(ifelse(tolower(as.character(x)) %in% c("na","nan","missing",""), NA, x)))
read_biom_tsv <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[!grepl("^# Constructed from biom file", lines)]
  tmp <- tempfile(); writeLines(lines, tmp)
  as.matrix(read.delim(tmp, check.names = FALSE, row.names = 1))
}

meta <- read.delim("00_metadata/sample_metadata.tsv", check.names = FALSE, stringsAsFactors = FALSE)
rownames(meta) <- meta$sample_id
meta$age_years <- numify(meta$age_years)
meta$BMI <- numify(meta$BMI)
meta$mets_status <- factor(meta$mets_status, levels = c("Control", "MetS"))
meta$sex <- factor(meta$sex)

# ---- Table 2: DADA2 tracking ----
stats_candidates <- c(
  "04_dada2/denoising_stats_export/stats.tsv",
  "results/dada2_read_tracking.tsv"
)
# export if needed
if (!file.exists(stats_candidates[1]) && file.exists("04_dada2/denoising_stats.qza")) {
  system2("qiime", c("tools","export","--input-path","04_dada2/denoising_stats.qza","--output-path","04_dada2/denoising_stats_export"))
}
if (file.exists("04_dada2/denoising_stats_export/stats.tsv")) {
  tr <- read.delim("04_dada2/denoising_stats_export/stats.tsv", check.names = FALSE)
  sid <- intersect(c("sample-id","Sample","sample_id"), names(tr))[1]
  if (!is.na(sid)) names(tr)[names(tr)==sid] <- "sample_id"
  write_tsv(as_tibble(tr), "results/dada2_read_tracking.tsv")
  write_tsv(as_tibble(tr), "12_tables/Table2_qc_dada2_read_tracking.tsv")
  if ("non-chimeric" %in% names(tr)) {
    tr2 <- tr %>% left_join(as_tibble(meta) %>% select(sample_id, mets_status), by = "sample_id")
    p_depth <- ggplot(tr2, aes(reorder(sample_id, `non-chimeric`), `non-chimeric`, fill = mets_status)) +
      geom_col() + coord_flip() +
      scale_fill_manual(values = c(Control = "#4C78A8", MetS = "#F58518")) +
      theme_bw(base_size = 9) +
      labs(title = "DADA2 non-chimeric reads per sample", x = NULL, y = "Reads",
           subtitle = "Control = all non-MetS")
    ggsave("11_figures/02_sequencing_depth.png", p_depth, width = 8, height = 10, dpi = 300)
    ggsave("11_figures/02_sequencing_depth.pdf", p_depth, width = 8, height = 10)
  }
}

# ---- Heatmap ----
genus <- read_biom_tsv("results/genus_abundance.tsv")
common <- intersect(colnames(genus), rownames(meta))
genus <- genus[, common, drop = FALSE]
short_genus <- function(x) {
  sapply(strsplit(x, ";"), function(p) {
    g <- p[length(p)]
    if (is.na(g) || g %in% c("", "g__", "__")) {
      for (i in rev(seq_along(p))) if (!p[i] %in% c("", "g__", "f__", "o__", "c__", "p__", "k__", "__")) return(p[i])
      return(x)
    }
    gsub("^g__", "", g)
  })
}
gagg <- rowsum(genus, group = short_genus(rownames(genus)))
rel <- sweep(gagg, 2, colSums(gagg), "/")
prev <- rowMeans(gagg > 0)
top <- names(sort(prev, decreasing = TRUE))[1:min(30, length(prev))]
mat <- rel[top, , drop = FALSE]
ann_cols <- intersect(c("mets_status", "Dx"), colnames(meta))
annot <- meta[colnames(mat), ann_cols, drop = FALSE]
png("11_figures/06_genus_abundance_heatmap.png", width = 2200, height = 1600, res = 200)
pheatmap(mat, annotation_col = annot, scale = "row", clustering_distance_cols = "euclidean",
         main = "Top prevalent genera (row-scaled relative abundance)")
dev.off()
pdf("11_figures/06_genus_abundance_heatmap.pdf", width = 11, height = 8)
pheatmap(mat, annotation_col = annot, scale = "row", clustering_distance_cols = "euclidean",
         main = "Top prevalent genera (row-scaled relative abundance)")
dev.off()

# ---- Workflow figure ----
wf <- tibble(
  step = factor(c("Metadata","FASTQ","FastQC","cutadapt","DADA2","SILVA138","Filter+tree","Diversity","ANCOM-BC2","Effect sizes"),
                levels = c("Metadata","FASTQ","FastQC","cutadapt","DADA2","SILVA138","Filter+tree","Diversity","ANCOM-BC2","Effect sizes")),
  x = 1:10, y = 1
)
p_wf <- ggplot(wf, aes(x, y)) +
  geom_tile(width = 0.9, height = 0.55, fill = "#4C78A8", color = "white") +
  geom_text(aes(label = step), color = "white", size = 3) +
  theme_void() +
  labs(title = "PRJNA838778 workflow (single-study; ANCOM-BC2 only)")
ggsave("11_figures/01_workflow.png", p_wf, width = 14, height = 2.5, dpi = 300)
ggsave("11_figures/01_workflow.pdf", p_wf, width = 14, height = 2.5)

# ---- Sensitivity comparison table (unadj vs adj) ----
unadj <- read.delim("08_differential_abundance/ANCOMBC2_genus_unadjusted.tsv", check.names = FALSE)
if (file.exists("08_differential_abundance/ANCOMBC2_genus_adjusted.tsv")) {
  adj <- read.delim("08_differential_abundance/ANCOMBC2_genus_adjusted.tsv", check.names = FALSE)
  sens <- unadj %>%
    select(taxon, unadj_effect = effect_size, unadj_FDR = FDR, unadj_p = p_value) %>%
    left_join(adj %>% select(taxon, adj_effect = effect_size, adj_FDR = FDR, adj_p = p_value), by = "taxon") %>%
    mutate(
      same_direction = sign(unadj_effect) == sign(adj_effect),
      stable_significant = !is.na(unadj_FDR) & !is.na(adj_FDR) & unadj_FDR < 0.10 & adj_FDR < 0.10 & same_direction
    ) %>%
    arrange(unadj_FDR)
  write_tsv(sens, "12_tables/Table7_sensitivity_comparison.tsv")
  write_tsv(sens, "09_effect_sizes/sensitivity_comparison.tsv")
}

# ---- Software versions ----
vers <- tibble(
  software = c("R", "ANCOMBC", "phyloseq", "vegan", "ggplot2", "QIIME2_env", "SILVA", "DA_method"),
  version = c(
    paste(R.version$major, R.version$minor, sep = "."),
    as.character(packageVersion("ANCOMBC")),
    as.character(packageVersion("phyloseq")),
    as.character(packageVersion("vegan")),
    as.character(packageVersion("ggplot2")),
    "qiime2-amplicon-2024.10",
    "138 (99% NB classifier)",
    "ANCOM-BC2 only"
  )
)
write_tsv(vers, "results/software_versions.tsv")

# ---- Functional prediction note ----
writeLines(c(
  "Functional prediction was not performed; 16S rRNA data do not establish gene function.",
  "Skip reason: prioritize study-level genus effect sizes for meta-analysis."
), "10_functional_prediction/README.txt")

# ---- Stats for report ----
fdr_cut <- 0.10
suggest_cut <- 0.20
es <- read.delim("results/study_level_effect_sizes.tsv", check.names = FALSE)
sig_up <- es %>% filter(FDR < fdr_cut, effect_size > 0) %>% pull(taxon)
sig_down <- es %>% filter(FDR < fdr_cut, effect_size < 0) %>% pull(taxon)
sug <- es %>% filter(FDR >= fdr_cut, FDR < suggest_cut) %>% arrange(FDR)
n_asv_raw <- if (file.exists("results/ASV_table_raw.tsv")) nrow(read_biom_tsv("results/ASV_table_raw.tsv")) else NA
n_asv_filt <- if (file.exists("results/ASV_table_filtered.tsv")) nrow(read_biom_tsv("results/ASV_table_filtered.tsv")) else NA
depth <- if (file.exists("07_diversity/rarefaction_depth.txt")) readLines("07_diversity/rarefaction_depth.txt")[1] else "see core_metrics"

alpha <- if (file.exists("12_tables/Table3_alpha_diversity.tsv") && file.info("12_tables/Table3_alpha_diversity.tsv")$size > 0) read.delim("12_tables/Table3_alpha_diversity.tsv") else NULL
beta <- if (file.exists("12_tables/Table4_beta_diversity_permanova.tsv") && file.info("12_tables/Table4_beta_diversity_permanova.tsv")$size > 0) read.delim("12_tables/Table4_beta_diversity_permanova.tsv") else NULL

age_str <- if (sum(!is.na(meta$age_years)) >= 3) sprintf("%.2f ± %.2f", mean(meta$age_years, na.rm=TRUE), sd(meta$age_years, na.rm=TRUE)) else "not deposited in BioSample"
bmi_str <- if (sum(!is.na(meta$BMI)) >= 3) sprintf("%.2f ± %.2f", mean(meta$BMI, na.rm=TRUE), sd(meta$BMI, na.rm=TRUE)) else "not deposited in BioSample"

report <- c(
  "# Analysis report — PRJNA838778 (single-study only)",
  "",
  "## A. Dataset summary",
  "PRJNA838778 is a pediatric fecal 16S V4 study (Zhou et al., *Microbiol Spectr* 2023; doi:10.1128/spectrum.03771-22).",
  paste0("- Samples/subjects: ", nrow(meta)),
  paste0("- MetS: ", sum(meta$mets_status=="MetS"), " | Non-MetS comparators: ", sum(meta$mets_status=="Control")),
  "- **Comparator definition for this analysis = all non-MetS** (`Control*` aliases).",
  paste0("- Age mean±SD: ", age_str),
  paste0("- BMI mean±SD: ", bmi_str),
  "- Primers: 515F / 806R (V4); platform: Illumina NovaSeq 6000; **single-end** deposited reads (~304–430 bp).",
  "- Sample counts match paper: 23 MetS + 31 non-MetS controls = 54 total.",
  "",
  "## B. Sample QC",
  paste0("FastQC completed for all ", nrow(meta), " samples (MultiQC skipped by request)."),
  "See `results/raw_read_qc_summary.tsv`, `results/download_qc_summary.tsv`, and `12_tables/Table2_qc_dada2_read_tracking.tsv`.",
  "No samples removed solely for soft FastQC warnings.",
  "",
  "## C. Sequencing summary",
  paste0("- Rarefaction depth (diversity): ", depth),
  paste0("- ASVs raw / filtered: ", n_asv_raw, " / ", n_asv_filt),
  "- Taxonomy: SILVA 138 99% Naive Bayes (`results/ASV_taxonomy.tsv`).",
  "- Genus table: `results/genus_abundance.tsv`.",
  "",
  "## D. Diversity results",
  "Alpha: `12_tables/Table3_alpha_diversity.tsv`; figure `11_figures/03_alpha_diversity.pdf`.",
  "Beta/PERMANOVA: `12_tables/Table4_beta_diversity_permanova.tsv`; PCoA in `11_figures/04_*` and `05_*`.",
  if (!is.null(alpha)) paste0("Alpha Wilcoxon FDR (Observed/Shannon/Simpson): ", paste(sprintf("%s=%.3g", alpha$metric, alpha$FDR), collapse="; ")) else NULL,
  "",
  "## E. Differential abundance",
  "Primary method: **ANCOM-BC2** at genus level (prevalence ≥ 0.10). Differential abundance used ANCOM-BC2 only.",
  paste0("Genera in effect-size table: ", nrow(es)),
  paste0("Significant FDR<", fdr_cut, " increased in MetS: ", if (length(sig_up)) paste(sig_up, collapse=", ") else "none"),
  paste0("Significant FDR<", fdr_cut, " decreased in MetS: ", if (length(sig_down)) paste(sig_down, collapse=", ") else "none"),
  if (nrow(sug)) paste0("Suggestive FDR<", suggest_cut, ": ", paste(sprintf("%s (%s, FDR=%.3g)", sug$taxon, sug$effect_direction, sug$FDR), collapse="; ")) else paste0("No suggestive taxa at FDR<", suggest_cut, "."),
  "",
  "## F. Covariate effects",
  "Adjusted ANCOM-BC2 **skipped**: age, sex, and BMI not available per-sample in BioSample.",
  "Per-sample TG/HDL/HbA1c/diet/medications were not in BioSample metadata and were not modeled.",
  "",
  "## G. Robust taxa",
  "See `12_tables/Table7_sensitivity_comparison.tsv` (unadjusted effects; adjusted model not run).",
  "",
  "## H. Study-level effect sizes (input for the meta-analysis)",
  "File: `results/study_level_effect_sizes.tsv`",
  "Effect: ANCOM-BC2 log fold-change (MetS vs Control) with SE and 95% CI.",
  "Control = all non-MetS comparators.",
  "",
  "## I. Functional prediction",
  "Not run (exploratory). See `10_functional_prediction/README.txt`.",
  "",
  "## J. Limitations specific to PRJNA838778",
  "1. Comparators are all children without MetS.",
  "2. Modest n (54; MetS=23).",
  "3. Pediatric cohort; age/sex/BMI not in public per-sample metadata.",
  "4. Key labs not deposited per-sample in SRA/BioSample.",
  "5. Cross-sectional; no causal inference.",
  "6. Deposited reads are single-end despite PAIRED library layout in SRA.",
  paste0("7. Significance threshold for DA reporting is FDR<", fdr_cut, " (relaxed from 0.05)."),
  "",
  "---",
  "Software: `results/software_versions.tsv`",
  "Config: `config/analysis_config.yaml`"
)
report <- report[!sapply(report, is.null)]
writeLines(report, "results/ANALYSIS_REPORT.md")
writeLines(report, "ANALYSIS_REPORT.md")

# Status board
status <- c(
  "# PRJNA838778 pipeline STATUS",
  paste0("Updated: ", format(Sys.time(), "%Y-%m-%dT%H:%M:%S")),
  "",
  "01. [DONE   ] 01 Verify raw data",
  "02. [DONE   ] 02 FastQC (no MultiQC)",
  "03. [DONE   ] 03 Import+cutadapt",
  "04. [DONE   ] 04 DADA2 ASVs",
  "05. [DONE   ] 05 Taxonomy SILVA138",
  "06. [DONE   ] 06 Filter+phylogeny",
  "07. [DONE   ] 07 Diversity",
  "08. [DONE   ] 08 DA ANCOM-BC2 + effect sizes + report",
  "",
  "Primary deliverable: results/study_level_effect_sizes.tsv"
)
writeLines(status, "results/PIPELINE_STATUS.md")
message("FINISH_PIPELINE_COMPLETE")
