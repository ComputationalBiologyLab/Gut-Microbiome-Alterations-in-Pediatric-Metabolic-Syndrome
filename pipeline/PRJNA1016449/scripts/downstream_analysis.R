#!/usr/bin/env Rscript
# PRJNA1016449 downstream analysis — single study only (no cross-study meta-analysis)
suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(ape)
  library(ggplot2)
  library(pheatmap)
  library(phyloseq)
})

root <- normalizePath(".")
dir.create("11_figures", showWarnings = FALSE)
dir.create("12_tables", showWarnings = FALSE)
dir.create("08_differential_abundance", showWarnings = FALSE)
dir.create("09_effect_sizes", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
dir.create("logs", showWarnings = FALSE)

bioproject <- "PRJNA1016449"
fdr_cut <- 0.10
suggest_cut <- 0.20
prev_cut <- 0.10

`%||%` <- function(a, b) if (!is.null(a) && length(a) && !is.na(a[1])) a else b

read_biom_tsv <- function(path) {
  # QIIME biom convert --to-tsv has a comment line then header with OTU ID
  lines <- readLines(path)
  lines <- lines[!grepl("^# Constructed from biom file", lines)]
  # First remaining line starts with #OTU ID or taxonomy
  tmp <- tempfile()
  writeLines(lines, tmp)
  df <- read.delim(tmp, check.names = FALSE, row.names = 1)
  as.matrix(df)
}

meta <- read.delim("00_metadata/sample_metadata.tsv", check.names = FALSE, stringsAsFactors = FALSE)
rownames(meta) <- meta$sample_id
# helper: metadata already has sample_id column — do not rownames_to_column again
meta_tbl <- function(m = meta) as_tibble(m)
numify <- function(x) suppressWarnings(as.numeric(ifelse(tolower(as.character(x)) %in% c("na","nan","missing",""), NA, x)))
meta$age_years <- numify(meta$age_years)
meta$BMI <- numify(meta$BMI)
meta$SBPI <- numify(meta$SBPI)
meta$DBPI <- numify(meta$DBPI)
meta$mets_status <- factor(meta$mets_status, levels = c("Control", "MetS"))
meta$clinical_group <- factor(meta$clinical_group, levels = c("Ob", "MetS"))
meta$BMIGroup <- factor(meta$BMIGroup, levels = c("BMI-A", "BMI-B"))
meta$sex <- factor(meta$sex)

# ---------- Table 1: dataset characteristics ----------
tab1 <- tibble(
  metric = c(
    "BioProject", "Publication", "Total_samples", "Total_subjects",
    "MetS_n", "Control_nonMetS_n",
    "Age_mean_sd", "Age_median_IQR",
    "Sex_male_n_pct", "Sex_female_n_pct",
    "BMI_mean_sd", "BMI_median_IQR",
    "BMIGroup_A_lt30", "BMIGroup_B_ge30",
    "Platform", "Region_16S", "Primers", "Read_layout", "Nominal_read_length",
    "Control_definition", "MetS_definition_summary"
  ),
  value = c(
    bioproject,
    "Gathungu et al. Obesities 2023; doi:10.3390/obesities3030021",
    nrow(meta), nrow(meta),
    sum(meta$mets_status == "MetS"), sum(meta$mets_status == "Control"),
    sprintf("%.2f ± %.2f", mean(meta$age_years), sd(meta$age_years)),
    sprintf("%.1f (%.1f–%.1f)", median(meta$age_years), quantile(meta$age_years, .25), quantile(meta$age_years, .75)),
    sprintf("%d (%.1f%%)", sum(meta$sex == "male"), 100 * mean(meta$sex == "male")),
    sprintf("%d (%.1f%%)", sum(meta$sex == "female"), 100 * mean(meta$sex == "female")),
    sprintf("%.2f ± %.2f", mean(meta$BMI), sd(meta$BMI)),
    sprintf("%.2f (%.2f–%.2f)", median(meta$BMI), quantile(meta$BMI, .25), quantile(meta$BMI, .75)),
    sum(meta$BMIGroup == "BMI-A"), sum(meta$BMIGroup == "BMI-B"),
    "Illumina MiSeq (600-cycle v3)", "V3–V4", "338F / 806R", "paired-end", "~300 bp",
    "Non-MetS comparator group (mets_status=Control)",
    "Obesity (>95th pct) + elevated SBP, TG, IFG, or hyperinsulinemia"
  )
)
write_tsv(tab1, "12_tables/Table1_dataset_characteristics.tsv")

# ---------- Helper: genus matrix ----------
genus_path <- "results/genus_abundance.tsv"
if (!file.exists(genus_path)) stop("Missing results/genus_abundance.tsv — run upstream steps first")
genus <- read_biom_tsv(genus_path)
# Align samples
common <- intersect(colnames(genus), rownames(meta))
genus <- genus[, common, drop = FALSE]
meta <- meta[common, , drop = FALSE]

# Clean genus names
rownames(genus) <- gsub(";__", ";", rownames(genus))
short_genus <- function(x) {
  sapply(strsplit(x, ";"), function(p) {
    g <- p[length(p)]
    if (is.na(g) || g == "" || g %in% c("g__", "__", "unclassified")) {
      # climb up
      for (i in rev(seq_along(p))) {
        if (!p[i] %in% c("", "g__", "f__", "o__", "c__", "p__", "k__", "__")) return(p[i])
      }
      return(x)
    }
    gsub("^g__", "", g)
  })
}
g_short <- short_genus(rownames(genus))
# Aggregate identical short names
genus_agg <- rowsum(genus, group = g_short)

rel <- sweep(genus_agg, 2, colSums(genus_agg), "/")
prev <- rowMeans(genus_agg > 0)
keep <- prev >= prev_cut
genus_da <- genus_agg[keep, , drop = FALSE]
rel_da <- rel[keep, , drop = FALSE]
message(sprintf("Genera for DA: %d / %d (prevalence >= %.2f)", sum(keep), nrow(genus_agg), prev_cut))

# ---------- Alpha diversity from filtered ASV table if available ----------
asv_path <- "results/ASV_table_filtered.tsv"
alpha <- NULL
if (file.exists(asv_path)) {
  asv <- read_biom_tsv(asv_path)
  asv <- asv[, common, drop = FALSE]
  # rarefy to depth file if present
  depth_file <- "07_diversity/rarefaction_depth.txt"
  if (file.exists(depth_file)) {
    depth <- as.integer(readLines(depth_file)[1])
  } else {
    depth <- floor(quantile(colSums(asv), 0.2))
  }
  keep_s <- colSums(asv) >= depth
  asv_r <- t(rrarefy(t(asv[, keep_s, drop = FALSE]), sample = depth))
  alpha <- tibble(
    sample_id = colnames(asv_r),
    observed = specnumber(t(asv_r)),
    shannon = diversity(t(asv_r), index = "shannon"),
    simpson = diversity(t(asv_r), index = "simpson")
  ) %>% left_join(meta_tbl() %>% select(sample_id, mets_status, clinical_group, age_years, sex, BMI, BMIGroup), by = "sample_id")
  write_tsv(alpha, "results/alpha_diversity_rarefied.tsv")
}

summarize_alpha <- function(df, metric) {
  x <- df[[metric]]
  g <- df$mets_status
  # Wilcoxon + Hodges-Lehmann-like median difference via quantile of pairwise diffs
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

tab3 <- NULL
if (!is.null(alpha)) {
  tab3 <- bind_rows(lapply(c("observed", "shannon", "simpson"), summarize_alpha, df = alpha))
  tab3$FDR <- p.adjust(tab3$p_value, method = "BH")
  write_tsv(tab3, "12_tables/Table3_alpha_diversity.tsv")

  # Figure: alpha diversity
  along <- alpha %>% pivot_longer(c(observed, shannon, simpson), names_to = "metric", values_to = "value")
  p_alpha <- ggplot(along, aes(mets_status, value, fill = mets_status)) +
    geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.8) +
    geom_jitter(width = 0.15, size = 1.5, alpha = 0.7) +
    facet_wrap(~metric, scales = "free_y") +
    scale_fill_manual(values = c(Control = "#4C78A8", MetS = "#F58518")) +
    theme_bw(base_size = 12) +
    labs(title = "Alpha diversity: MetS vs non-MetS (PRJNA1016449)",
         subtitle = sprintf("Rarefaction depth = %s; Control = non-MetS comparators", if (exists("depth")) depth else "NA"),
         x = NULL, y = NULL) +
    theme(legend.position = "none")
  ggsave("11_figures/03_alpha_diversity.png", p_alpha, width = 9, height = 4.5, dpi = 300)
  ggsave("11_figures/03_alpha_diversity.pdf", p_alpha, width = 9, height = 4.5)
}

# ---------- Beta diversity ----------
read_dist <- function(path) {
  as.dist(as.matrix(read.delim(path, row.names = 1, check.names = FALSE)))
}

beta_rows <- list()
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
    labs(title = title,
         x = sprintf("PCoA1 (%.1f%%)", var1),
         y = sprintf("PCoA2 (%.1f%%)", var2),
         color = "Group")
  ggsave(paste0(outfile, ".png"), p, width = 6.5, height = 5, dpi = 300)
  ggsave(paste0(outfile, ".pdf"), p, width = 6.5, height = 5)
  invisible(pc)
}

run_permanova <- function(dist_obj, meta_df, metric_name) {
  md <- meta_df[labels(dist_obj), , drop = FALSE]
  # Unadjusted
  set.seed(1016449)
  fit1 <- adonis2(dist_obj ~ mets_status, data = md, permutations = 999)
  # Dispersion
  bd <- betadisper(dist_obj, md$mets_status)
  set.seed(1016449)
  bd_p <- permutest(bd, permutations = 999)
  # Adjusted if variables available / not confounded
  # Check collinearity: mets_status vs BMIGroup association
  fit2 <- tryCatch({
    set.seed(1016449)
    adonis2(dist_obj ~ mets_status + age_years + sex + BMI, data = md, permutations = 999, by = "margin")
  }, error = function(e) NULL)

  extract_term <- function(fit, term, model_label) {
    if (is.null(fit)) return(NULL)
    rn <- rownames(fit)
    i <- which(rn == term)
    if (!length(i)) return(NULL)
    tibble(
      metric = metric_name,
      model = model_label,
      term = term,
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
  bind_rows(
    extract_term(fit1, "mets_status", "mets_status"),
    extract_term(fit2, "mets_status", "mets_status + age_years + sex + BMI")
  )
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
  if (!file.exists(pth)) next
  d <- read_dist(pth)
  # align
  labs <- intersect(labels(d), rownames(meta))
  m <- as.matrix(d)[labs, labs]
  d2 <- as.dist(m)
  md2 <- meta[labs, , drop = FALSE]
  tab4[[nm]] <- run_permanova(d2, md2, nm)
  if (nm == "braycurtis") plot_pcoa(d2, md2, "Bray–Curtis PCoA (PRJNA1016449)", "11_figures/04_pcoa_braycurtis")
  if (nm == "weighted_unifrac") plot_pcoa(d2, md2, "Weighted UniFrac PCoA (PRJNA1016449)", "11_figures/05_pcoa_weighted_unifrac")
  if (nm == "unweighted_unifrac") plot_pcoa(d2, md2, "Unweighted UniFrac PCoA (PRJNA1016449)", "11_figures/05b_pcoa_unweighted_unifrac")
}
if (length(tab4)) {
  tab4 <- bind_rows(tab4)
  write_tsv(tab4, "12_tables/Table4_beta_diversity_permanova.tsv")
}

# ---------- Sequencing depth figure ----------
if (file.exists("results/dada2_read_tracking.tsv")) {
  tr <- read.delim("results/dada2_read_tracking.tsv", check.names = FALSE)
  # QIIME export may have sample-id column
  sid_col <- intersect(c("sample-id", "Sample", "sample_id"), names(tr))[1]
  if (!is.na(sid_col)) {
    colnames(tr)[colnames(tr) == sid_col] <- "sample_id"
    depth_col <- intersect(c("non-chimeric", "non_chimeric", "input"), names(tr))
    # plot non-chimeric if present
    if ("non-chimeric" %in% names(tr)) {
      tr2 <- tr %>% left_join(meta_tbl() %>% select(sample_id, mets_status), by = "sample_id")
      p_depth <- ggplot(tr2, aes(reorder(sample_id, `non-chimeric`), `non-chimeric`, fill = mets_status)) +
        geom_col() + coord_flip() +
        scale_fill_manual(values = c(Control = "#4C78A8", MetS = "#F58518")) +
        theme_bw(base_size = 9) +
        labs(title = "DADA2 non-chimeric reads per sample", x = NULL, y = "Reads")
      ggsave("11_figures/02_sequencing_depth.png", p_depth, width = 8, height = 10, dpi = 300)
      ggsave("11_figures/02_sequencing_depth.pdf", p_depth, width = 8, height = 10)
    }
  }
  write_tsv(as_tibble(tr), "12_tables/Table2_qc_dada2_read_tracking.tsv")
}

# ---------- ANCOM-BC2 ----------
ancom_ok <- requireNamespace("ANCOMBC", quietly = TRUE)
maaslin_ok <- FALSE  # single DA method only: ANCOM-BC2

tax_tab <- data.frame(Genus = rownames(genus_da), row.names = rownames(genus_da))
otu <- otu_table(genus_da, taxa_are_rows = TRUE)
sam <- sample_data(meta)
TAX <- tax_table(as.matrix(tax_tab))
ps <- phyloseq(otu, sam, TAX)

ancom_res <- NULL
ancom_adj <- NULL
if (ancom_ok) {
  message("Running ANCOM-BC2 (unadjusted)...")
  set.seed(1016449)
  out1 <- ANCOMBC::ancombc2(
    data = ps,
    assay_name = "counts",
    tax_level = NULL,
    fix_formula = "mets_status",
    rand_formula = NULL,
    p_adj_method = "BH",
    prv_cut = 0,  # already filtered
    lib_cut = 0,
    s0_perc = 0.05,
    group = "mets_status",
    struc_zero = TRUE,
    neg_lb = TRUE,
    alpha = 0.05,
    n_cl = 1,
    verbose = TRUE
  )
  res1 <- out1$res
  # Column names vary by ANCOMBC version; detect MetS effect columns
  lfc_col <- grep("lfc.*MetS|MetS.*lfc", names(res1), value = TRUE, ignore.case = TRUE)
  se_col <- grep("se.*MetS|MetS.*se", names(res1), value = TRUE, ignore.case = TRUE)
  p_col <- grep("^p_.*MetS|MetS.*[^q]p$|p_mets_statusMetS", names(res1), value = TRUE, ignore.case = TRUE)
  q_col <- grep("^q_.*MetS|MetS.*q$|q_mets_statusMetS", names(res1), value = TRUE, ignore.case = TRUE)
  # Fallback: columns containing mets_statusMetS
  if (!length(lfc_col)) lfc_col <- grep("mets_statusMetS", names(res1), value = TRUE)[1]
  if (!length(se_col)) se_col <- grep("^se_", names(res1), value = TRUE)[1]
  if (!length(p_col)) p_col <- grep("^p_", names(res1), value = TRUE)[1]
  if (!length(q_col)) q_col <- grep("^q_", names(res1), value = TRUE)[1]
  message("ANCOMBC columns used: ", paste(c(lfc_col[1], se_col[1], p_col[1], q_col[1]), collapse=", "))

  ancom_res <- res1 %>%
    transmute(
      taxon = taxon %||% rownames(res1),
      taxonomic_level = "genus",
      effect_size = .data[[lfc_col[1]]],
      standard_error = .data[[se_col[1]]],
      p_value = .data[[p_col[1]]],
      FDR = .data[[q_col[1]]]
    ) %>%
    mutate(
      lower_CI = effect_size - 1.96 * standard_error,
      upper_CI = effect_size + 1.96 * standard_error,
      effect_direction = case_when(
        is.na(effect_size) ~ "NA",
        effect_size > 0 ~ "increased_in_MetS",
        effect_size < 0 ~ "decreased_in_MetS",
        TRUE ~ "null"
      ),
      model = "unadjusted",
      covariates = "none",
      analysis_method = "ANCOM-BC2",
      prevalence = prev[taxon],
      MetS_n = sum(meta$mets_status == "MetS"),
      Control_n = sum(meta$mets_status == "Control"),
      mean_rel_MetS = rowMeans(rel_da[taxon, meta$mets_status == "MetS", drop = FALSE]),
      mean_rel_Control = rowMeans(rel_da[taxon, meta$mets_status == "Control", drop = FALSE]),
      median_rel_MetS = apply(rel_da[taxon, meta$mets_status == "MetS", drop = FALSE], 1, median),
      median_rel_Control = apply(rel_da[taxon, meta$mets_status == "Control", drop = FALSE], 1, median)
    )

  # Adjusted model
  message("Running ANCOM-BC2 (adjusted: age + sex + BMI)...")
  set.seed(1016449)
  out2 <- tryCatch(ANCOMBC::ancombc2(
    data = ps,
    assay_name = "counts",
    tax_level = NULL,
    fix_formula = "mets_status + age_years + sex + BMI",
    rand_formula = NULL,
    p_adj_method = "BH",
    prv_cut = 0,
    lib_cut = 0,
    group = "mets_status",
    struc_zero = TRUE,
    neg_lb = TRUE,
    alpha = 0.05,
    n_cl = 1,
    verbose = TRUE
  ), error = function(e) { message("Adjusted ANCOM-BC2 failed: ", e$message); NULL })

  if (!is.null(out2)) {
    res2 <- out2$res
    lfc_col <- grep("lfc_mets_statusMetS|lfc.*MetS", names(res2), value = TRUE, ignore.case = TRUE)[1]
    se_col <- grep("se_mets_statusMetS|^se_.*MetS", names(res2), value = TRUE, ignore.case = TRUE)[1]
    p_col <- grep("^p_mets_statusMetS|^p_.*MetS", names(res2), value = TRUE, ignore.case = TRUE)[1]
    q_col <- grep("^q_mets_statusMetS|^q_.*MetS", names(res2), value = TRUE, ignore.case = TRUE)[1]
    ancom_adj <- res2 %>%
      transmute(
        taxon = taxon %||% rownames(res2),
        taxonomic_level = "genus",
        effect_size = .data[[lfc_col]],
        standard_error = .data[[se_col]],
        p_value = .data[[p_col]],
        FDR = .data[[q_col]],
        lower_CI = effect_size - 1.96 * .data[[se_col]],
        upper_CI = effect_size + 1.96 * .data[[se_col]],
        effect_direction = ifelse(effect_size > 0, "increased_in_MetS", "decreased_in_MetS"),
        model = "adjusted",
        covariates = "age_years,sex,BMI",
        analysis_method = "ANCOM-BC2",
        prevalence = prev[taxon],
        MetS_n = sum(meta$mets_status == "MetS"),
        Control_n = sum(meta$mets_status == "Control")
      )
  }
  write_tsv(ancom_res, "08_differential_abundance/ANCOMBC2_genus_unadjusted.tsv")
  if (!is.null(ancom_adj)) write_tsv(ancom_adj, "08_differential_abundance/ANCOMBC2_genus_adjusted.tsv")
} else {
  message("ANCOMBC not installed — skipping primary DA")
}

# ---------- MaAsLin2 secondary ----------
maas_res <- NULL
if (maaslin_ok) {
  message("Running MaAsLin2...")
  md_m <- meta_tbl()
  # relative abundance input; TSS already
  fit_dir <- "08_differential_abundance/maaslin2_unadjusted"
  unlink(fit_dir, recursive = TRUE)
  fit <- Maaslin2::Maaslin2(
    input_data = as.data.frame(t(rel_da)),
    input_metadata = md_m %>% column_to_rownames("sample_id"),
    output = fit_dir,
    fixed_effects = c("mets_status"),
    reference = c("mets_status,Control"),
    normalization = "NONE",
    transform = "LOG",
    analysis_method = "LM",
    correction = "BH",
    standardize = FALSE,
    min_prevalence = 0,
    plot_heatmap = FALSE,
    plot_scatter = FALSE
  )
  maas_res <- read.delim(file.path(fit_dir, "all_results.tsv"), check.names = FALSE) %>%
    filter(metadata == "mets_status") %>%
    transmute(
      taxon = feature,
      taxonomic_level = "genus",
      effect_size = coef,
      standard_error = stderr,
      p_value = pval,
      FDR = qval,
      lower_CI = coef - 1.96 * stderr,
      upper_CI = coef + 1.96 * stderr,
      effect_direction = ifelse(coef > 0, "increased_in_MetS", "decreased_in_MetS"),
      model = "unadjusted",
      covariates = "none",
      analysis_method = "MaAsLin2_LM_log",
      prevalence = prev[feature],
      MetS_n = sum(meta$mets_status == "MetS"),
      Control_n = sum(meta$mets_status == "Control")
    )
  write_tsv(maas_res, "08_differential_abundance/Maaslin2_genus_unadjusted.tsv")
}

# ---------- Primary effect-size table for meta-analysis ----------
if (!is.null(ancom_res)) {
  es <- ancom_res %>%
    mutate(BioProject = bioproject) %>%
    select(BioProject, taxon, taxonomic_level, effect_direction, effect_size, standard_error,
           lower_CI, upper_CI, p_value, FDR, prevalence, MetS_n, Control_n, model, covariates, analysis_method,
           mean_rel_MetS, mean_rel_Control, median_rel_MetS, median_rel_Control)
  write_tsv(es, "results/study_level_effect_sizes.tsv")
  write_tsv(es, "09_effect_sizes/study_level_effect_sizes.tsv")
  write_tsv(es, "12_tables/Table6_primary_study_level_effect_sizes.tsv")

  # Classification table
  class_tbl <- es %>%
    mutate(
      classification = case_when(
        is.na(FDR) ~ "insufficient_data",
        FDR < fdr_cut & effect_size > 0 ~ "significantly_increased_in_MetS",
        FDR < fdr_cut & effect_size < 0 ~ "significantly_decreased_in_MetS",
        FDR < suggest_cut ~ "suggestive_nonsignificant",
        TRUE ~ "non_significant"
      )
    ) %>%
    select(taxon, direction = effect_direction, effect = effect_size, SE = standard_error, FDR, classification, p_value)
  # Add insufficient prevalence genera
  insuff <- tibble(
    taxon = setdiff(names(prev), rownames(genus_da)),
    direction = NA_character_,
    effect = NA_real_,
    SE = NA_real_,
    FDR = NA_real_,
    classification = "insufficient_prevalence",
    p_value = NA_real_
  )
  class_tbl <- bind_rows(class_tbl, insuff)
  write_tsv(class_tbl, "09_effect_sizes/genus_classification.tsv")
  write_tsv(class_tbl, "results/genus_classification.tsv")

  tab5 <- es %>% arrange(FDR)
  write_tsv(tab5, "12_tables/Table5_differential_abundance_ANCOMBC2.tsv")

  # Volcano
  vol <- es %>% mutate(sig = FDR < fdr_cut)
  p_vol <- ggplot(vol, aes(effect_size, -log10(p_value), color = sig)) +
    geom_point(alpha = 0.85, size = 2) +
    scale_color_manual(values = c(`FALSE` = "grey70", `TRUE` = "#E45756")) +
    theme_bw(base_size = 12) +
    labs(title = "Genus differential abundance (ANCOM-BC2)",
         subtitle = "MetS vs non-MetS | PRJNA1016449",
         x = "log fold-change (MetS / Control)", y = expression(-log[10](p)), color = "FDR<0.10")
  if (any(vol$sig, na.rm = TRUE) && requireNamespace("ggrepel", quietly = TRUE)) {
    p_vol <- p_vol + ggrepel::geom_text_repel(data = subset(vol, sig), aes(label = taxon), size = 3, max.overlaps = 20)
  }
  ggsave("11_figures/07_volcano_ancombc2.png", p_vol, width = 7, height = 5.5, dpi = 300)
  ggsave("11_figures/07_volcano_ancombc2.pdf", p_vol, width = 7, height = 5.5)

  # Top differential genera bar
  top <- vol %>% arrange(FDR) %>% head(20)
  p_top <- ggplot(top, aes(reorder(taxon, effect_size), effect_size, fill = effect_direction)) +
    geom_col() + coord_flip() +
    scale_fill_manual(values = c(increased_in_MetS = "#E45756", decreased_in_MetS = "#4C78A8")) +
    theme_bw(base_size = 11) +
    labs(title = "Top genera by ANCOM-BC2 FDR", x = NULL, y = "logFC (MetS / Control)", fill = NULL)
  ggsave("11_figures/08_top_differential_genera.png", p_top, width = 8, height = 6, dpi = 300)
  ggsave("11_figures/08_top_differential_genera.pdf", p_top, width = 8, height = 6)
}

# ---------- Heatmap top genera ----------
top_prev <- names(sort(prev, decreasing = TRUE))[1:min(30, length(prev))]
mat <- rel[top_prev, , drop = FALSE]
annot <- meta %>% select(mets_status, BMIGroup, sex)
png("11_figures/06_genus_abundance_heatmap.png", width = 2200, height = 1600, res = 200)
pheatmap(mat, annotation_col = annot, scale = "row", clustering_distance_cols = "euclidean",
         main = "Top prevalent genera (relative abundance, row-scaled)")
dev.off()
pdf("11_figures/06_genus_abundance_heatmap.pdf", width = 11, height = 8)
pheatmap(mat, annotation_col = annot, scale = "row", clustering_distance_cols = "euclidean",
         main = "Top prevalent genera (relative abundance, row-scaled)")
dev.off()

# ---------- Sensitivity: prevalence thresholds ----------
sens_list <- list()
if (ancom_ok) {
  for (pc in c(0.05, 0.10, 0.20)) {
    g2 <- genus_agg[prev >= pc, , drop = FALSE]
    if (nrow(g2) < 5) next
    ps2 <- phyloseq(otu_table(g2, taxa_are_rows = TRUE), sample_data(meta),
                    tax_table(as.matrix(data.frame(Genus = rownames(g2), row.names = rownames(g2)))))
    o <- tryCatch(ANCOMBC::ancombc2(
      data = ps2, fix_formula = "mets_status", p_adj_method = "BH", prv_cut = 0,
      group = "mets_status", struc_zero = TRUE, neg_lb = TRUE, alpha = 0.05, n_cl = 1, verbose = FALSE
    ), error = function(e) NULL)
    if (is.null(o)) next
    rr <- o$res
    lfc_col <- grep("lfc_mets_statusMetS|lfc.*MetS", names(rr), value = TRUE, ignore.case = TRUE)[1]
    q_col <- grep("^q_mets_statusMetS|^q_.*MetS", names(rr), value = TRUE, ignore.case = TRUE)[1]
    sens_list[[as.character(pc)]] <- rr %>%
      transmute(taxon = taxon, effect_size = .data[[lfc_col]], FDR = .data[[q_col]], prevalence_threshold = pc)
  }
}
sens <- if (length(sens_list)) bind_rows(sens_list) else NULL
if (!is.null(sens) && !is.null(ancom_res)) {
  wide <- sens %>%
    mutate(sig = FDR < fdr_cut) %>%
    select(taxon, prevalence_threshold, sig, effect_size) %>%
    pivot_wider(names_from = prevalence_threshold, values_from = c(sig, effect_size), names_glue = "{.value}_prev{prevalence_threshold}")
  # Stability vs adjusted / MaAsLin2
  stab <- ancom_res %>% select(taxon, unadj_effect = effect_size, unadj_FDR = FDR)
  if (!is.null(ancom_adj)) stab <- stab %>% left_join(ancom_adj %>% select(taxon, adj_effect = effect_size, adj_FDR = FDR), by = "taxon")
  if (!is.null(maas_res)) stab <- stab %>% left_join(maas_res %>% select(taxon, maas_effect = effect_size, maas_FDR = FDR), by = "taxon")
  stab <- stab %>% left_join(wide, by = "taxon") %>%
    mutate(stable_sig_unadj_adj = if ("adj_FDR" %in% names(.)) (unadj_FDR < fdr_cut & adj_FDR < fdr_cut & sign(unadj_effect) == sign(adj_effect)) else NA)
  write_tsv(stab, "12_tables/Table7_sensitivity_comparison.tsv")
  write_tsv(stab, "09_effect_sizes/sensitivity_comparison.tsv")
}

# ---------- Family-level exploratory ----------
fam_path <- "results/family_abundance.tsv"
if (file.exists(fam_path) && ancom_ok) {
  fam <- read_biom_tsv(fam_path)
  fam <- fam[, common, drop = FALSE]
  fshort <- short_genus(rownames(fam)) # reuse last rank parser
  fam_agg <- rowsum(fam, group = fshort)
  fprev <- rowMeans(fam_agg > 0)
  fam_da <- fam_agg[fprev >= prev_cut, , drop = FALSE]
  ps_f <- phyloseq(otu_table(fam_da, taxa_are_rows = TRUE), sample_data(meta),
                   tax_table(as.matrix(data.frame(Family = rownames(fam_da), row.names = rownames(fam_da)))))
  of <- tryCatch(ANCOMBC::ancombc2(
    data = ps_f, fix_formula = "mets_status", p_adj_method = "BH", prv_cut = 0,
    group = "mets_status", struc_zero = TRUE, neg_lb = TRUE, alpha = 0.05, n_cl = 1, verbose = FALSE
  ), error = function(e) NULL)
  if (!is.null(of)) {
    rr <- of$res
    lfc_col <- grep("lfc_mets_statusMetS|lfc.*MetS", names(rr), value = TRUE, ignore.case = TRUE)[1]
    se_col <- grep("se_mets_statusMetS|^se_.*MetS", names(rr), value = TRUE, ignore.case = TRUE)[1]
    p_col <- grep("^p_mets_statusMetS|^p_.*MetS", names(rr), value = TRUE, ignore.case = TRUE)[1]
    q_col <- grep("^q_mets_statusMetS|^q_.*MetS", names(rr), value = TRUE, ignore.case = TRUE)[1]
    fam_out <- rr %>% transmute(
      BioProject = bioproject, taxon = taxon, taxonomic_level = "family",
      effect_size = .data[[lfc_col]], standard_error = .data[[se_col]],
      lower_CI = effect_size - 1.96 * standard_error, upper_CI = effect_size + 1.96 * standard_error,
      p_value = .data[[p_col]], FDR = .data[[q_col]],
      effect_direction = ifelse(effect_size > 0, "increased_in_MetS", "decreased_in_MetS"),
      prevalence = fprev[taxon], MetS_n = sum(meta$mets_status == "MetS"), Control_n = sum(meta$mets_status == "Control"),
      model = "unadjusted", covariates = "none", analysis_method = "ANCOM-BC2"
    )
    write_tsv(fam_out, "08_differential_abundance/ANCOMBC2_family_unadjusted.tsv")
  }
}

# ---------- Workflow diagram (simple ggplot) ----------
wf <- tibble(
  step = factor(c("Metadata\nSRA/BioSample", "FASTQ\nacquisition", "FastQC/\nMultiQC", "cutadapt\nprimer trim",
                  "DADA2\nASVs", "SILVA138\ntaxonomy", "Filter +\nphylogeny", "Diversity", "ANCOM-BC2\n+ MaAsLin2",
                  "Study-level\neffect sizes"), levels = c("Metadata\nSRA/BioSample", "FASTQ\nacquisition", "FastQC/\nMultiQC", "cutadapt\nprimer trim",
                  "DADA2\nASVs", "SILVA138\ntaxonomy", "Filter +\nphylogeny", "Diversity", "ANCOM-BC2\n+ MaAsLin2",
                  "Study-level\neffect sizes")),
  x = 1:10, y = 1
)
p_wf <- ggplot(wf, aes(x, y)) +
  geom_tile(width = 0.9, height = 0.5, fill = "#4C78A8", color = "white") +
  geom_text(aes(label = step), color = "white", size = 2.7, lineheight = 0.9) +
  geom_segment(data = tibble(x = 1:9 + 0.45, xend = 2:10 - 0.45, y = 1, yend = 1),
               aes(x = x, xend = xend, y = y, yend = yend),
               inherit.aes = FALSE, arrow = arrow(length = unit(0.15, "cm")), linewidth = 0.4) +
  theme_void() +
  labs(title = "PRJNA1016449 analysis workflow (single-study; no meta-analysis yet)")
ggsave("11_figures/01_workflow.png", p_wf, width = 14, height = 2.8, dpi = 300)
ggsave("11_figures/01_workflow.pdf", p_wf, width = 14, height = 2.8)

# ---------- Software versions ----------
vers <- tibble(
  software = c("R", "vegan", "phyloseq", "ggplot2",
               "ANCOMBC", "Maaslin2", "QIIME2_env"),
  version = c(
    paste(R.version$major, R.version$minor, sep = "."),
    as.character(packageVersion("vegan")),
    as.character(packageVersion("phyloseq")),
    as.character(packageVersion("ggplot2")),
    if (ancom_ok) as.character(packageVersion("ANCOMBC")) else "not_installed",
    if (maaslin_ok) as.character(packageVersion("Maaslin2")) else "not_installed",
    "qiime2-amplicon-2024.10"
  )
)
write_tsv(vers, "results/software_versions.tsv")

# ---------- Analysis report ----------
sig_up <- if (!is.null(ancom_res)) ancom_res %>% filter(FDR < fdr_cut, effect_size > 0) %>% pull(taxon) else character()
sig_down <- if (!is.null(ancom_res)) ancom_res %>% filter(FDR < fdr_cut, effect_size < 0) %>% pull(taxon) else character()
robust <- if (exists("stab") && is.data.frame(stab) && "stable_sig_unadj_adj" %in% names(stab)) {
  stab %>% filter(stable_sig_unadj_adj) %>% pull(taxon)
} else character()

report <- c(
  "# Analysis report — PRJNA1016449 (single-study only)",
  "",
  "## A. Dataset summary",
  "PRJNA1016449 is a pediatric fecal 16S V3–V4 MiSeq study (Gathungu et al., Obesities 2023) comparing children with obesity **with metabolic syndrome (MetS)** to children with obesity **without MetS (Ob)**.",
  paste0("- Samples/subjects: ", nrow(meta)),
  paste0("- MetS: ", sum(meta$mets_status=="MetS"), " | Non-MetS comparators: ", sum(meta$mets_status=="Control")),
  "- **Comparator group for all MetS contrasts = non-MetS** (`mets_status = Control`).",
  paste0("- Age mean±SD: ", sprintf("%.2f ± %.2f", mean(meta$age_years), sd(meta$age_years))),
  paste0("- BMI mean±SD: ", sprintf("%.2f ± %.2f", mean(meta$BMI), sd(meta$BMI))),
  "- Primers: 338F / 806R (V3–V4); platform: Illumina MiSeq 2×300.",
  "- BioSample Dx counts (MetS=16, Ob=34) differ slightly from paper Table 1 (17/33); BioSample labels used.",
  "",
  "## B. Sample QC",
  "See `results/raw_read_qc_summary.tsv` and `results/dada2_read_tracking.tsv`.",
  "No samples were removed solely for soft FastQC warnings; DADA2 retention is sample-wise in Table 2.",
  "",
  "## C. Sequencing summary",
  "ASV features after filtering are in `results/ASV_table_filtered.tsv`.",
  "Taxonomy: SILVA 138 99% Naive Bayes (`results/ASV_taxonomy.tsv`).",
  "Genus table: `results/genus_abundance.tsv`.",
  "",
  "## D. Diversity results",
  "Alpha: `12_tables/Table3_alpha_diversity.tsv` and `11_figures/03_alpha_diversity.pdf`.",
  "Beta/PERMANOVA/PERMDISP: `12_tables/Table4_beta_diversity_permanova.tsv` and PCoA figures.",
  "",
  "## E. Differential abundance",
  paste0("Primary method: **ANCOM-BC2** at genus level (prevalence ≥ ", prev_cut, ")."),
  "Secondary DA method: skipped (single-method analysis).",
  paste0("Significant genera (FDR < ", fdr_cut, "):"),
  paste0("- Increased in MetS: ", if (length(sig_up)) paste(sig_up, collapse=", ") else "none"),
  paste0("- Decreased in MetS: ", if (length(sig_down)) paste(sig_down, collapse=", ") else "none"),
  "",
  "## F. Covariate effects",
  "Adjusted ANCOM-BC2 formula: mets_status + age_years + sex + BMI (when model converged).",
  "Compare `ANCOMBC2_genus_unadjusted.tsv` vs `ANCOMBC2_genus_adjusted.tsv`.",
  "Diet/medication/T2DM per-sample covariates were **not deposited** in BioSample and were not modeled.",
  "",
  "## G. Robust taxa",
  paste0("Taxa significant in both unadjusted and adjusted models (same direction): ", if (length(robust)) paste(robust, collapse=", ") else "see Table7 / none"),
  "",
  "## H. Study-level effect sizes (input for the meta-analysis)",
  "File: `results/study_level_effect_sizes.tsv`",
  "Effect definition: ANCOM-BC2 log fold-change (MetS vs Control) with SE and 95% CI.",
  "Control = non-MetS comparator group.",
  "Ready for random-effects meta-analysis with other BioProjects **after those studies are processed separately**.",
  "",
  "## I. Functional prediction",
  "Exploratory only; run PICRUSt2 separately if `picrust2` env is used. 16S-inferred functions are not causal proof of gene content.",
  "",
  "## J. Limitations specific to PRJNA1016449",
  "1. Comparators are children without MetS; interpret pooled effects as MetS versus non-MetS.",
  "2. Modest sample size (n=50; MetS≈16).",
  "3. Pediatric cohort (ages ~10–19); may not generalize to adult MetS meta-analysis without stratification.",
  "4. Key labs (TG, HDL, HbA1c) reported in paper but not in public per-sample metadata.",
  "5. Cross-sectional; no causal inference.",
  "6. BioSample Dx counts disagree slightly with published Table 1.",
  "7. Barcoded primer constructs required cutadapt removal of upstream barcode+primer.",
  "",
  "---",
  "Software versions: `results/software_versions.tsv`",
  "Config: `config/analysis_config.yaml`"
)
writeLines(report, "results/ANALYSIS_REPORT.md")
writeLines(report, "ANALYSIS_REPORT.md")
message("Wrote ANALYSIS_REPORT.md")
message("Downstream analysis complete.")
