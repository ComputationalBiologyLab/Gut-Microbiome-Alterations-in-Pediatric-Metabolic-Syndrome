# Sensitivity of alpha- and beta-diversity results to rarefaction depth .
# Re-runs the diversity tests at (a) the depth used in the paper and (b) the depth of the shallowest sample, which retains every sample.
suppressPackageStartupMessages(library(vegan))
B <- paste0(sub("/$", "", Sys.getenv("MS_ROOT", ".")), "/")
coh <- list(US=list(id="PRJNA1016449", seed=1016449, adj=TRUE), China=list(id="PRJNA838778", seed=838778, adj=FALSE), Mexico=list(id="PRJNA819279", seed=819279, adj=FALSE))
res <- list()
for (nm in names(coh)) {
  cc <- coh[[nm]]; d <- paste0(B, cc$id, "_analysis/")
  asv <- read.delim(paste0(d, "results/ASV_table_filtered.tsv"), skip = 1, row.names = 1, check.names = FALSE)
  md <- read.delim(paste0(d, "00_metadata/sample_metadata.tsv"), check.names = FALSE)
  md <- md[!duplicated(md$sample_id), ]; rownames(md) <- md$sample_id
  md <- md[colnames(asv), ]
  md$grp <- factor(ifelse(md$mets_status == "MetS", "MetS", "Control"), levels = c("Control", "MetS"))
  if (cc$adj) { md$age_years <- as.numeric(md$age_years); md$BMI <- as.numeric(md$BMI) }
  depth_used <- as.integer(readLines(paste0(d, "07_diversity/rarefaction_depth.txt"))[1])
  for (dl in c("depth used in the paper", "shallowest sample (all samples retained)")) {
    depth <- if (dl == "depth used in the paper") depth_used else min(colSums(asv))
    keep <- colSums(asv) >= depth; set.seed(cc$seed)
    R <- t(rrarefy(t(asv[, keep, drop = FALSE]), sample = depth)); g <- md$grp[keep]
    a <- data.frame(observed = specnumber(t(R)), shannon = diversity(t(R), "shannon"), simpson = diversity(t(R), "simpson"))
    pw <- sapply(names(a), function(m) wilcox.test(a[[m]] ~ g, exact = FALSE)$p.value); fdr <- p.adjust(pw, "BH")
    for (m in names(a)) res[[length(res) + 1]] <- data.frame(Cohort = nm, Depth = dl, Rarefaction_depth = depth, Samples_kept = sum(keep), MetS_kept = sum(g == "MetS"), Comparator_kept = sum(g == "Control"),
      Test = paste("alpha:", m), Estimate_MetS_minus_comparator = median(a[[m]][g == "MetS"]) - median(a[[m]][g == "Control"]), R2 = NA, P = pw[[m]], FDR_BH = fdr[[m]], Betadisper_P = NA)
    for (dist in c("bray", "jaccard")) {
      D <- vegdist(t(R), method = dist, binary = (dist == "jaccard")); mdk <- md[keep, ]; mdk$g <- g
      set.seed(cc$seed)
      f <- if (cc$adj) adonis2(D ~ g + age_years + sex + BMI, data = mdk, permutations = 999, by = "margin") else adonis2(D ~ g, data = mdk, permutations = 999)
      i <- which(rownames(f) %in% c("g", "Model"))[1]; bd <- permutest(betadisper(D, g), permutations = 999)
      res[[length(res) + 1]] <- data.frame(Cohort = nm, Depth = dl, Rarefaction_depth = depth, Samples_kept = sum(keep), MetS_kept = sum(g == "MetS"), Comparator_kept = sum(g == "Control"),
        Test = paste("PERMANOVA:", dist), Estimate_MetS_minus_comparator = NA, R2 = f$R2[i], P = f$`Pr(>F)`[i], FDR_BH = NA, Betadisper_P = bd$tab$`Pr(>F)`[1])
    }
  }
}
out <- do.call(rbind, res)
write.table(out, "rarefaction_sensitivity.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
print(out, digits = 3)
