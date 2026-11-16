#!/usr/bin/env Rscript
# Phylum-level ANCOM-BC2 for all three BioProjects (same method as genus pipeline)
suppressPackageStartupMessages({
  library(tidyverse)
  library(phyloseq)
  library(ANCOMBC)
})

run_study <- function(study_dir, bioproject, prev_cut = 0.10) {
  setwd(study_dir)
  message("\n=== ", bioproject, " ===")

  read_biom_tsv <- function(path) {
    lines <- readLines(path, warn = FALSE)
    lines <- lines[!grepl("^# Constructed from biom file", lines)]
    tmp <- tempfile()
    writeLines(lines, tmp)
    as.matrix(read.delim(tmp, check.names = FALSE, row.names = 1))
  }

  extract_phylum <- function(x) {
    sapply(strsplit(x, ";"), function(p) {
      hit <- grep("^p__", p, value = TRUE)
      if (length(hit)) return(sub("^p__", "", hit[1]))
      "Unassigned"
    })
  }

  meta <- read.delim("00_metadata/sample_metadata.tsv", check.names = FALSE, stringsAsFactors = FALSE)
  rownames(meta) <- meta$sample_id
  meta$mets_status <- factor(meta$mets_status, levels = c("Control", "MetS"))

  genus <- read_biom_tsv("results/genus_abundance.tsv")
  common <- intersect(colnames(genus), rownames(meta))
  genus <- genus[, common, drop = FALSE]
  meta <- meta[common, , drop = FALSE]

  phyla <- extract_phylum(rownames(genus))
  phylum_agg <- rowsum(genus, group = phyla)
  rel <- sweep(phylum_agg, 2, colSums(phylum_agg), "/")
  prev <- rowMeans(phylum_agg > 0)
  keep <- prev >= prev_cut
  phylum_da <- phylum_agg[keep, , drop = FALSE]
  rel_da <- rel[keep, , drop = FALSE]
  message(sprintf("Phyla for DA: %d / %d (prevalence >= %.2f)", sum(keep), nrow(phylum_agg), prev_cut))

  ps <- phyloseq(
    otu_table(phylum_da, taxa_are_rows = TRUE),
    sample_data(meta),
    tax_table(as.matrix(data.frame(Phylum = rownames(phylum_da), row.names = rownames(phylum_da))))
  )

  set.seed(as.integer(gsub("^PRJNA", "", bioproject)))
  out <- ancombc2(
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
    verbose = FALSE,
    pseudo_sens = FALSE
  )
  res <- out$res
  pick_col <- function(df, patterns) {
    for (pat in patterns) {
      hit <- grep(pat, names(df), value = TRUE, ignore.case = TRUE)
      if (length(hit)) return(hit[1])
    }
    NA_character_
  }
  lfc_col <- pick_col(res, c("^lfc_mets_statusMetS$", "lfc_.*MetS"))
  se_col <- pick_col(res, c("^se_mets_statusMetS$", "se_.*MetS"))
  p_col <- pick_col(res, c("^p_mets_statusMetS$", "^p_.*MetS"))
  q_col <- pick_col(res, c("^q_mets_statusMetS$", "^q_.*MetS"))
  tax_col <- if ("taxon" %in% names(res)) "taxon" else names(res)[1]

  es <- res %>%
    transmute(
      BioProject = bioproject,
      taxon = .data[[tax_col]],
      taxonomic_level = "phylum",
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
      Control_n = sum(meta$mets_status == "Control")
    ) %>%
    arrange(FDR, p_value)

  out_path <- "results/study_level_effect_sizes_phylum.tsv"
  write_tsv(es, out_path)
  message("Wrote ", out_path, " (n=", nrow(es), ")")
  sig <- es %>% filter(FDR < 0.1)
  if (nrow(sig)) print(sig %>% select(taxon, effect_size, p_value, FDR))
  es
}

# MS_ROOT: directory that contains the three <BioProject>_analysis folders (default: current directory)
root <- Sys.getenv("MS_ROOT", getwd())
studies <- list(
  list(dir = file.path(root, "PRJNA1016449_analysis"), bp = "PRJNA1016449"),
  list(dir = file.path(root, "PRJNA838778_analysis"), bp = "PRJNA838778"),
  list(dir = file.path(root, "PRJNA819279_analysis"), bp = "PRJNA819279")
)

all_es <- lapply(studies, function(s) run_study(s$dir, s$bp))
combined <- bind_rows(all_es)
dir.create(file.path(root, "meta_analysis"), showWarnings = FALSE)
write_tsv(combined, file.path(root, "meta_analysis", "phylum_study_level_effects_all.tsv"))
message("\nCombined: ", nrow(combined), " rows")
