#!/usr/bin/env Rscript
# PRJNA838778: alpha/beta diversity tables + figures (mirrors PRJNA1016449 deliverables)
suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(ape)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("--file=", "", args[grep("--file=", args)])
if (length(file_arg)) {
  root <- normalizePath(file.path(dirname(file_arg), ".."))
  setwd(root)
}

bioproject <- "PRJNA838778"
dir.create("11_figures", showWarnings = FALSE)
dir.create("12_tables", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

`%||%` <- function(a, b) if (!is.null(a) && length(a) && !is.na(a[1])) a else b

read_biom_tsv <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[!grepl("^# Constructed from biom file", lines)]
  tmp <- tempfile()
  writeLines(lines, tmp)
  as.matrix(read.delim(tmp, check.names = FALSE, row.names = 1))
}

meta <- read.delim("00_metadata/sample_metadata.tsv", check.names = FALSE, stringsAsFactors = FALSE)
rownames(meta) <- meta$sample_id
meta_tbl <- function(m = meta) as_tibble(m)
numify <- function(x) suppressWarnings(as.numeric(ifelse(tolower(as.character(x)) %in% c("na", "nan", "missing", ""), NA, x)))
meta$age_years <- numify(meta$age_years)
meta$BMI <- numify(meta$BMI)
meta$mets_status <- factor(meta$mets_status, levels = c("Control", "MetS"))
meta$sex <- factor(meta$sex)

age_mean <- mean(meta$age_years, na.rm = TRUE)
age_sd <- sd(meta$age_years, na.rm = TRUE)
bmi_mean <- mean(meta$BMI, na.rm = TRUE)
bmi_sd <- sd(meta$BMI, na.rm = TRUE)

tab1 <- tibble(
  metric = c(
    "BioProject", "Publication", "Total_samples", "Total_subjects",
    "MetS_n", "Control_all_non_MetS_n",
    "Age_mean_sd", "Age_median_IQR",
    "Sex_male_n_pct", "Sex_female_n_pct",
    "BMI_mean_sd", "BMI_median_IQR",
    "Platform", "Region_16S", "Primers", "Read_layout", "Nominal_read_length",
    "Control_definition", "MetS_definition_summary"
  ),
  value = c(
    bioproject,
    "Zhou et al. Microbiol Spectr 2023; doi:10.1128/spectrum.03771-22",
    nrow(meta), nrow(meta),
    sum(meta$mets_status == "MetS"), sum(meta$mets_status == "Control"),
    if (is.finite(age_mean)) sprintf("%.2f ± %.2f", age_mean, age_sd) else "not in BioSample",
    if (sum(!is.na(meta$age_years)) >= 3) sprintf("%.1f (%.1f–%.1f)", median(meta$age_years, na.rm = TRUE), quantile(meta$age_years, .25, na.rm = TRUE), quantile(meta$age_years, .75, na.rm = TRUE)) else "not in BioSample",
    if (sum(meta$sex == "male", na.rm = TRUE)) sprintf("%d (%.1f%%)", sum(meta$sex == "male", na.rm = TRUE), 100 * mean(meta$sex == "male", na.rm = TRUE)) else "not in BioSample",
    if (sum(meta$sex == "female", na.rm = TRUE)) sprintf("%d (%.1f%%)", sum(meta$sex == "female", na.rm = TRUE), 100 * mean(meta$sex == "female", na.rm = TRUE)) else "not in BioSample",
    if (is.finite(bmi_mean)) sprintf("%.2f ± %.2f", bmi_mean, bmi_sd) else "not in BioSample",
    if (sum(!is.na(meta$BMI)) >= 3) sprintf("%.2f (%.2f–%.2f)", median(meta$BMI, na.rm = TRUE), quantile(meta$BMI, .25, na.rm = TRUE), quantile(meta$BMI, .75, na.rm = TRUE)) else "not in BioSample",
    "Illumina NovaSeq 6000 (Novogene)", "V4", "515F / 806R", "single-end", "~304–430 bp (primers removed in deposit)",
    "All non-MetS samples (Control*) used as controls (mets_status=Control)",
    "Pediatric MetS per Chinese Medical Society pediatric guidelines (Zhou et al. 2023)",
    "54 samples: 23 MetS (MS*) + 31 Control (Control*); matches paper"
  )
)
write_tsv(tab1, "12_tables/Table1_dataset_characteristics.tsv")

# download_qc_summary.tsv (single-end analogue of PRJNA1016449)
manifest <- read.delim("00_metadata/fastq_manifest.tsv", check.names = FALSE, stringsAsFactors = FALSE)
dl <- meta %>%
  select(sample_id, run_accession = run_accession, BioSample = biosample) %>%
  left_join(manifest %>% rename(sample_id = `sample-id`), by = "sample_id") %>%
  mutate(
    Run = run_accession,
    r1_exists = file.exists(`absolute-filepath`),
    r2_exists = FALSE,
    r1_bytes = ifelse(r1_exists, file.info(`absolute-filepath`)$size, NA_real_),
    r2_bytes = NA_real_,
    paired = FALSE,
    duplicate_run = duplicated(run_accession) | duplicated(run_accession, fromLast = TRUE)
  ) %>%
  select(sample_id, Run, BioSample, r1_exists, r2_exists, r1_bytes, r2_bytes, paired, duplicate_run)
write_tsv(dl, "results/download_qc_summary.tsv")

asv_path <- "results/ASV_table_filtered.tsv"
stopifnot(file.exists(asv_path))

common <- intersect(colnames(read_biom_tsv(asv_path)), rownames(meta))
meta <- meta[common, , drop = FALSE]
asv <- read_biom_tsv(asv_path)
asv <- asv[, common, drop = FALSE]

depth_file <- "07_diversity/rarefaction_depth.txt"
depth <- if (file.exists(depth_file)) as.integer(readLines(depth_file)[1]) else floor(quantile(colSums(asv), 0.2))
keep_s <- colSums(asv) >= depth
asv_r <- t(rrarefy(t(asv[, keep_s, drop = FALSE]), sample = depth))
alpha <- tibble(
  sample_id = colnames(asv_r),
  observed = specnumber(t(asv_r)),
  shannon = diversity(t(asv_r), index = "shannon"),
  simpson = diversity(t(asv_r), index = "simpson")
) %>%
  left_join(meta_tbl() %>% select(sample_id, mets_status, age_years, sex, BMI), by = "sample_id")
write_tsv(alpha, "results/alpha_diversity_rarefied.tsv")

summarize_alpha <- function(df, metric) {
  x <- df[[metric]]
  g <- df$mets_status
  w <- wilcox.test(x ~ g, conf.int = TRUE, exact = FALSE)
  splitx <- split(x, g)
  tibble(
    metric = metric,
    test = "Wilcoxon rank-sum",
    Control_n = length(splitx$Control),
    MetS_n = length(splitx$MetS),
    Control_median = median(splitx$Control),
    Control_IQR = IQR(splitx$Control),
    MetS_median = median(splitx$MetS),
    MetS_IQR = IQR(splitx$MetS),
    effect_estimate = unname(w$estimate %||% (median(splitx$MetS) - median(splitx$Control))),
    effect_name = "pseudo_median_diff_MetS_minus_Control",
    CI_low = w$conf.int[1],
    CI_high = w$conf.int[2],
    p_value = w$p.value
  )
}

tab3 <- bind_rows(lapply(c("observed", "shannon", "simpson"), summarize_alpha, df = alpha))
tab3$FDR <- p.adjust(tab3$p_value, method = "BH")
write_tsv(tab3, "12_tables/Table3_alpha_diversity.tsv")

along <- alpha %>% pivot_longer(c(observed, shannon, simpson), names_to = "metric", values_to = "value")
p_alpha <- ggplot(along, aes(mets_status, value, fill = mets_status)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.7) +
  facet_wrap(~metric, scales = "free_y") +
  scale_fill_manual(values = c(Control = "#4C78A8", MetS = "#F58518")) +
  theme_bw(base_size = 12) +
  labs(
    title = "Alpha diversity: MetS vs Control (PRJNA838778)",
    subtitle = sprintf("Rarefaction depth = %s; Control = all non-MetS (n=%d kept)", depth, sum(keep_s)),
    x = NULL, y = NULL
  ) +
  theme(legend.position = "none")
ggsave("11_figures/03_alpha_diversity.png", p_alpha, width = 9, height = 4.5, dpi = 300)
ggsave("11_figures/03_alpha_diversity.pdf", p_alpha, width = 9, height = 4.5)

read_dist <- function(path) as.dist(as.matrix(read.delim(path, row.names = 1, check.names = FALSE)))

plot_pcoa <- function(dist_obj, meta_df, title, outfile) {
  pc <- pcoa(dist_obj)
  vec <- pc$vectors[, 1:2, drop = FALSE]
  df <- data.frame(sample_id = rownames(vec), PC1 = vec[, 1], PC2 = vec[, 2]) %>%
    left_join(meta_tbl(meta_df), by = "sample_id")
  var1 <- 100 * pc$values$Relative_eig[1]
  var2 <- 100 * pc$values$Relative_eig[2]
  p <- ggplot(df, aes(PC1, PC2, color = mets_status)) +
    geom_point(size = 2.8, alpha = 0.9) +
    stat_ellipse(level = 0.95, linewidth = 0.4) +
    scale_color_manual(values = c(Control = "#4C78A8", MetS = "#F58518")) +
    theme_bw(base_size = 12) +
    labs(
      title = title,
      x = sprintf("PCoA1 (%.1f%%)", var1),
      y = sprintf("PCoA2 (%.1f%%)", var2),
      color = "Group"
    )
  ggsave(paste0(outfile, ".png"), p, width = 6.5, height = 5, dpi = 300)
  ggsave(paste0(outfile, ".pdf"), p, width = 6.5, height = 5)
}

run_permanova <- function(dist_obj, meta_df, metric_name) {
  ord <- labels(dist_obj)
  md <- meta_df[match(ord, rownames(meta_df)), , drop = FALSE]
  rownames(md) <- ord
  md$mets_status <- droplevels(factor(md$mets_status, levels = c("Control", "MetS")))
  set.seed(838778)
  fit1 <- adonis2(dist_obj ~ mets_status, data = md, permutations = 999)
  bd <- betadisper(dist_obj, md$mets_status)
  set.seed(838778)
  bd_p <- permutest(bd, permutations = 999)
  fit2 <- NULL

  extract_term <- function(fit, model_label) {
    if (is.null(fit)) return(NULL)
    rn <- rownames(fit)
    i <- grep("mets_status", rn, ignore.case = TRUE)[1]
    term_label <- if (!is.na(i)) rn[i] else "mets_status"
    if (is.na(i)) i <- which(rn == "Model")[1]
    if (is.na(i)) i <- which(!rn %in% c("Total", "Residual"))[1]
    if (is.na(i)) return(NULL)
    tibble(
      metric = metric_name,
      model = model_label,
      term = term_label,
      R2 = fit$R2[i],
      F_stat = fit$F[i],
      df = fit$Df[i],
      p_value = fit$`Pr(>F)`[i],
      permutations = 999,
      n_samples = nrow(md),
      MetS_n = sum(md$mets_status == "MetS"),
      Control_n = sum(md$mets_status == "Control"),
      betadisper_F = bd_p$tab$F[1],
      betadisper_p = bd_p$tab$`Pr(>F)`[1],
      package = "vegan::adonis2 / betadisper"
    )
  }
  out <- extract_term(fit1, "mets_status")
  if (is.null(out) || nrow(out) == 0) {
    message("adonis2 rownames (", metric_name, "): ", paste(rownames(fit1), collapse = ", "))
  }
  out
}

dist_map <- list(
  braycurtis = "07_diversity/bray_curtis_export/distance-matrix.tsv",
  jaccard = "07_diversity/jaccard_export/distance-matrix.tsv",
  unweighted_unifrac = "07_diversity/unweighted_unifrac_export/distance-matrix.tsv",
  weighted_unifrac = "07_diversity/weighted_unifrac_export/distance-matrix.tsv"
)

tab4 <- list()
for (nm in names(dist_map)) {
  pth <- dist_map[[nm]]
  if (!file.exists(pth)) {
    message("Missing distance matrix: ", pth)
    next
  }
  message("PERMANOVA: ", nm, " (", pth, ")")
  d <- tryCatch(read_dist(pth), error = function(e) { message("read_dist failed: ", e$message); NULL })
  if (is.null(d)) next
  labs <- intersect(labels(d), rownames(meta))
  if (length(labs) < 4) {
    message("Too few overlapping samples for ", nm, ": ", length(labs))
    next
  }
  d2 <- as.dist(as.matrix(d)[labs, labs])
  md2 <- meta[labs, , drop = FALSE]
  tab4[[nm]] <- tryCatch(run_permanova(d2, md2, nm), error = function(e) {
    message("PERMANOVA failed for ", nm, ": ", e$message)
    NULL
  })
  if (nm == "braycurtis") try(plot_pcoa(d2, md2, "Bray–Curtis PCoA (PRJNA838778)", "11_figures/04_pcoa_braycurtis"), silent = TRUE)
  if (nm == "weighted_unifrac") try(plot_pcoa(d2, md2, "Weighted UniFrac PCoA (PRJNA838778)", "11_figures/05_pcoa_weighted_unifrac"), silent = TRUE)
  if (nm == "unweighted_unifrac") try(plot_pcoa(d2, md2, "Unweighted UniFrac PCoA (PRJNA838778)", "11_figures/05b_pcoa_unweighted_unifrac"), silent = TRUE)
}
tab4_df <- bind_rows(tab4)
if (nrow(tab4_df)) {
  write_tsv(tab4_df, "12_tables/Table4_beta_diversity_permanova.tsv")
} else {
  message("WARNING: no PERMANOVA rows written — check distance-matrix exports")
}

# Table7: unadjusted only (adjusted skipped — no covariates in BioSample)
if (file.exists("08_differential_abundance/ANCOMBC2_genus_unadjusted.tsv")) {
  unadj <- read.delim("08_differential_abundance/ANCOMBC2_genus_unadjusted.tsv", check.names = FALSE)
  tab7 <- unadj %>%
    transmute(
      taxon,
      unadj_effect = effect_size,
      unadj_FDR = FDR,
      unadj_p = p_value,
      adj_effect = NA_real_,
      adj_FDR = NA_real_,
      adj_p = NA_real_,
      same_direction = NA,
      stable_significant = FALSE,
      note = "adjusted_ANCOMBC2_skipped_no_covariates_in_BioSample"
    ) %>%
    arrange(unadj_FDR)
  write_tsv(tab7, "12_tables/Table7_sensitivity_comparison.tsv")
  write_tsv(tab7, "09_effect_sizes/sensitivity_comparison.tsv")
}

message("DIVERSITY_EXPORTS_COMPLETE")
