#!/usr/bin/env Rscript
# Fast single-method DA: ANCOM-BC2 only (MetS vs Control=non-MetS (T2DM excluded))
suppressPackageStartupMessages({
  library(tidyverse)
  library(phyloseq)
  library(ANCOMBC)
})

dir.create("08_differential_abundance", showWarnings = FALSE)
dir.create("09_effect_sizes", showWarnings = FALSE)
dir.create("12_tables", showWarnings = FALSE)
dir.create("11_figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

bioproject <- "PRJNA819279"
fdr_cut <- 0.10
suggest_cut <- 0.20
prev_cut <- 0.10

read_biom_tsv <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[!grepl("^# Constructed from biom file", lines)]
  tmp <- tempfile(); writeLines(lines, tmp)
  as.matrix(read.delim(tmp, check.names = FALSE, row.names = 1))
}

numify <- function(x) suppressWarnings(as.numeric(ifelse(tolower(as.character(x)) %in% c("na","nan","missing",""), NA, x)))

meta <- read.delim("00_metadata/sample_metadata.tsv", check.names = FALSE, stringsAsFactors = FALSE)
rownames(meta) <- meta$sample_id
meta$age_years <- numify(meta$age_years)
meta$BMI <- numify(meta$BMI)
meta$mets_status <- factor(meta$mets_status, levels = c("Control", "MetS"))
meta$sex <- factor(meta$sex)
stopifnot(all(as.character(meta$mets_status) %in% c("Control", "MetS")))
stopifnot(!any(grepl("T2DM", as.character(meta$mets_status), ignore.case = TRUE)))

genus <- read_biom_tsv("results/genus_abundance.tsv")
common <- intersect(colnames(genus), rownames(meta))
genus <- genus[, common, drop = FALSE]
meta <- meta[common, , drop = FALSE]

short_genus <- function(x) {
  sapply(strsplit(x, ";"), function(p) {
    g <- p[length(p)]
    if (is.na(g) || g %in% c("", "g__", "__", "unclassified")) {
      for (i in rev(seq_along(p))) {
        if (!p[i] %in% c("", "g__", "f__", "o__", "c__", "p__", "k__", "__")) return(p[i])
      }
      return(x)
    }
    gsub("^g__", "", g)
  })
}
genus_agg <- rowsum(genus, group = short_genus(rownames(genus)))
rel <- sweep(genus_agg, 2, colSums(genus_agg), "/")
prev <- rowMeans(genus_agg > 0)
keep <- prev >= prev_cut
genus_da <- genus_agg[keep, , drop = FALSE]
rel_da <- rel[keep, , drop = FALSE]
message(sprintf("Genera for DA: %d / %d (prevalence >= %.2f)", sum(keep), nrow(genus_agg), prev_cut))

# ANCOM-BC2 fails on zero-variance taxa (overall or within any group)
row_var <- apply(genus_da, 1, var)
zero_var <- rownames(genus_da)[!is.finite(row_var) | row_var <= 0]
for (grp in levels(meta$mets_status)) {
  sub <- genus_da[, meta$mets_status == grp, drop = FALSE]
  gv <- apply(sub, 1, var)
  zero_var <- unique(c(zero_var, rownames(sub)[!is.finite(gv) | gv <= 0]))
}
if (length(zero_var)) {
  message("Removing zero-variance genera: ", paste(zero_var, collapse = ", "))
  genus_da <- genus_da[!rownames(genus_da) %in% zero_var, , drop = FALSE]
  rel_da <- rel_da[rownames(genus_da), , drop = FALSE]
}
message(sprintf("Genera after zero-variance filter: %d", nrow(genus_da)))

ps <- phyloseq(
  otu_table(genus_da, taxa_are_rows = TRUE),
  sample_data(meta),
  tax_table(as.matrix(data.frame(Genus = rownames(genus_da), row.names = rownames(genus_da))))
)

message("Running ANCOM-BC2 (unadjusted; pseudo_sens = FALSE)...")
set.seed(819279)
out1 <- ancombc2(
  data = ps,
  assay_name = "counts",
  tax_level = NULL,
  fix_formula = "mets_status",
  rand_formula = NULL,
  p_adj_method = "BH",
  prv_cut = 0,
  lib_cut = 0,
  s0_perc = 0.05,
  group = "mets_status",
  struc_zero = TRUE,
  neg_lb = TRUE,
  alpha = 0.05,
  n_cl = 1,
  verbose = TRUE,
  pseudo_sens = FALSE
)
res1 <- out1$res
message("ANCOMBC columns: ", paste(names(res1), collapse = ", "))

pick_col <- function(df, patterns) {
  for (pat in patterns) {
    hit <- grep(pat, names(df), value = TRUE, ignore.case = TRUE)
    if (length(hit)) return(hit[1])
  }
  NA_character_
}
lfc_col <- pick_col(res1, c("^lfc_mets_statusMetS$", "lfc_.*MetS", "lfc_mets_status"))
se_col  <- pick_col(res1, c("^se_mets_statusMetS$", "se_.*MetS", "^se_mets_status"))
p_col   <- pick_col(res1, c("^p_mets_statusMetS$", "^p_.*MetS", "^p_mets_status"))
q_col   <- pick_col(res1, c("^q_mets_statusMetS$", "^q_.*MetS", "^q_mets_status"))
tax_col <- if ("taxon" %in% names(res1)) "taxon" else names(res1)[1]
message("Using cols: ", paste(c(tax_col, lfc_col, se_col, p_col, q_col), collapse = " | "))
stopifnot(!is.na(lfc_col), !is.na(se_col), !is.na(p_col), !is.na(q_col))

es <- res1 %>%
  transmute(
    BioProject = bioproject,
    taxon = .data[[tax_col]],
    taxonomic_level = "genus",
    effect_size = .data[[lfc_col]],
    standard_error = .data[[se_col]],
    p_value = .data[[p_col]],
    FDR = .data[[q_col]],
    lower_CI = effect_size - 1.96 * standard_error,
    upper_CI = effect_size + 1.96 * standard_error,
    effect_direction = case_when(
      is.na(effect_size) ~ NA_character_,
      effect_size > 0 ~ "increased_in_MetS",
      effect_size < 0 ~ "decreased_in_MetS",
      TRUE ~ "null"
    ),
    model = "unadjusted",
    covariates = "none",
    analysis_method = "ANCOM-BC2",
    prevalence = as.numeric(prev[as.character(taxon)]),
    MetS_n = sum(meta$mets_status == "MetS"),
    Control_n = sum(meta$mets_status == "Control"),
    mean_rel_MetS = as.numeric(rowMeans(rel_da[as.character(taxon), meta$mets_status == "MetS", drop = FALSE])),
    mean_rel_Control = as.numeric(rowMeans(rel_da[as.character(taxon), meta$mets_status == "Control", drop = FALSE])),
    median_rel_MetS = as.numeric(apply(rel_da[as.character(taxon), meta$mets_status == "MetS", drop = FALSE], 1, median)),
    median_rel_Control = as.numeric(apply(rel_da[as.character(taxon), meta$mets_status == "Control", drop = FALSE], 1, median))
  ) %>%
  arrange(FDR, p_value)

write_tsv(es, "08_differential_abundance/ANCOMBC2_genus_unadjusted.tsv")
write_tsv(es, "results/study_level_effect_sizes.tsv")
write_tsv(es, "09_effect_sizes/study_level_effect_sizes.tsv")
write_tsv(es, "12_tables/Table5_differential_abundance_ANCOMBC2.tsv")
write_tsv(es, "12_tables/Table6_primary_study_level_effect_sizes.tsv")

class_tbl <- es %>%
  mutate(classification = case_when(
    is.na(FDR) ~ "insufficient_data",
    FDR < fdr_cut & effect_size > 0 ~ "significantly_increased_in_MetS",
    FDR < fdr_cut & effect_size < 0 ~ "significantly_decreased_in_MetS",
    FDR < suggest_cut ~ "suggestive_nonsignificant",
    TRUE ~ "non_significant"
  )) %>%
  select(taxon, direction = effect_direction, effect = effect_size, SE = standard_error, FDR, classification, p_value)
insuff <- tibble(
  taxon = setdiff(names(prev), rownames(genus_da)),
  direction = NA_character_, effect = NA_real_, SE = NA_real_, FDR = NA_real_,
  classification = "insufficient_prevalence", p_value = NA_real_
)
write_tsv(bind_rows(class_tbl, insuff), "09_effect_sizes/genus_classification.tsv")
write_tsv(bind_rows(class_tbl, insuff), "results/genus_classification.tsv")

# Adjusted model skipped: age/sex/BMI not deposited per-sample for this BioProject
message("Skipping adjusted ANCOM-BC2 (no age/sex/BMI in BioSample metadata).")
out2 <- NULL
writeLines(
  "Adjusted ANCOM-BC2 skipped: age_years, sex, and BMI not available per-sample in BioSample.",
  "08_differential_abundance/ANCOMBC2_adjusted_SKIPPED.txt"
)

if (!is.null(out2)) {
  res2 <- out2$res
  lfc2 <- pick_col(res2, c("^lfc_mets_statusMetS$", "lfc_.*MetS"))
  se2 <- pick_col(res2, c("^se_mets_statusMetS$", "se_.*MetS"))
  p2 <- pick_col(res2, c("^p_mets_statusMetS$", "^p_.*MetS"))
  q2 <- pick_col(res2, c("^q_mets_statusMetS$", "^q_.*MetS"))
  tax2 <- if ("taxon" %in% names(res2)) "taxon" else names(res2)[1]
  adj <- res2 %>% transmute(
    taxon = .data[[tax2]], effect_size = .data[[lfc2]], standard_error = .data[[se2]],
    p_value = .data[[p2]], FDR = .data[[q2]],
    lower_CI = effect_size - 1.96 * standard_error,
    upper_CI = effect_size + 1.96 * standard_error,
    effect_direction = ifelse(effect_size > 0, "increased_in_MetS", "decreased_in_MetS"),
    model = "adjusted", covariates = "age_years,sex,BMI", analysis_method = "ANCOM-BC2"
  )
  write_tsv(adj, "08_differential_abundance/ANCOMBC2_genus_adjusted.tsv")
}

# Simple volcano + top plot (label FDR-significant taxa)
has_repel <- requireNamespace("ggrepel", quietly = TRUE)
if (has_repel) suppressPackageStartupMessages(library(ggrepel))
vol <- es %>% mutate(
  sig = !is.na(FDR) & FDR < fdr_cut,
  label = ifelse(sig, as.character(taxon), NA_character_),
  y = -log10(pmax(p_value, 1e-300))
)
p_vol <- ggplot(vol, aes(effect_size, y, color = sig)) +
  geom_point(alpha = 0.85, size = 2.2) +
  scale_color_manual(values = c(`FALSE` = "grey70", `TRUE` = "#E45756"),
                     labels = c(`FALSE` = paste0("FDR ≥ ", fdr_cut), `TRUE` = paste0("FDR < ", fdr_cut))) +
  {
    if (has_repel) {
      geom_text_repel(aes(label = label), size = 3.2, color = "black",
                      max.overlaps = Inf, box.padding = 0.4, point.padding = 0.25,
                      min.segment.length = 0, seed = 1, na.rm = TRUE, show.legend = FALSE)
    } else {
      geom_text(aes(label = label), size = 3.2, color = "black",
                vjust = -0.7, na.rm = TRUE, show.legend = FALSE)
    }
  } +
  theme_bw(base_size = 12) +
  labs(title = "Genus DA — ANCOM-BC2 (PRJNA819279)",
       subtitle = paste0("MetS vs Control (non-MetS; T2DM excluded); labeled if FDR < ", fdr_cut),
       x = "logFC (MetS / Control)", y = expression(-log[10](p)), color = NULL)
ggsave("11_figures/07_volcano_ancombc2.png", p_vol, width = 7.5, height = 5.5, dpi = 300)
ggsave("11_figures/07_volcano_ancombc2.pdf", p_vol, width = 7.5, height = 5.5)

top <- vol %>% arrange(FDR) %>% head(20)
p_top <- ggplot(top, aes(reorder(taxon, effect_size), effect_size, fill = effect_direction)) +
  geom_col() + coord_flip() +
  scale_fill_manual(values = c(increased_in_MetS = "#E45756", decreased_in_MetS = "#4C78A8", `NA` = "grey60"), na.value = "grey60") +
  theme_bw(base_size = 11) +
  labs(title = "Top genera by ANCOM-BC2 FDR", x = NULL, y = "logFC (MetS / Control)", fill = NULL)
ggsave("11_figures/08_top_differential_genera.png", p_top, width = 8, height = 6, dpi = 300)
ggsave("11_figures/08_top_differential_genera.pdf", p_top, width = 8, height = 6)

sig_up <- es %>% filter(FDR < fdr_cut, effect_size > 0) %>% pull(taxon)
sig_down <- es %>% filter(FDR < fdr_cut, effect_size < 0) %>% pull(taxon)
message("Significant increased in MetS: ", if (length(sig_up)) paste(sig_up, collapse = ", ") else "none")
message("Significant decreased in MetS: ", if (length(sig_down)) paste(sig_down, collapse = ", ") else "none")
message("Wrote results/study_level_effect_sizes.tsv (n=", nrow(es), ")")
message("DA COMPLETE")
