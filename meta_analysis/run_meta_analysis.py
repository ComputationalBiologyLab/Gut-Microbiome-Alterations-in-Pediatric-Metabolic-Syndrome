#!/usr/bin/env python3
"""Three-cohort random-effects meta-analysis of genus-level ANCOM-BC2 effects.

Inputs (one table per cohort, produced by the single-study pipeline):
    <MS_ROOT>/<BioProject>_analysis/results/study_level_effect_sizes.tsv

Steps:
    1. Keep genera that are present in all three cohorts with a prevalence of at least 10% and a finite
       effect size and standard error. With three cohorts a majority direction always exists, so the
       "interpretable" criterion in `select_taxa_three` never excludes a genus.
    2. DerSimonian-Laird random-effects pooling, Benjamini-Hochberg FDR, Cochran Q, tau2 and I2.
    3. Reproducibility classification (strongly concordant, directionally concordant, heterogeneous,
       no clear association).

Outputs are written to <MS_ROOT>/meta_analysis/.

Usage:
    MS_ROOT=/path/to/project python meta_analysis/run_meta_analysis.py
"""

from __future__ import annotations

import math
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(os.environ.get("MS_ROOT", ".")).resolve()
META = ROOT / "meta_analysis"
REAL_101 = ROOT / "PRJNA1016449_analysis" / "results" / "study_level_effect_sizes.tsv"
REAL_838 = ROOT / "PRJNA838778_analysis" / "results" / "study_level_effect_sizes.tsv"
REAL_819 = ROOT / "PRJNA819279_analysis" / "results" / "study_level_effect_sizes.tsv"

MIN_PREVALENCE = 0.10
Z95 = 1.959963984540054


def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    adj = np.empty(n, dtype=float)
    prev = 1.0
    for i in range(n - 1, -1, -1):
        rank = i + 1
        val = ranked[i] * n / rank
        prev = min(prev, val)
        adj[i] = prev
    out = np.empty(n, dtype=float)
    out[order] = np.clip(adj, 0, 1)
    return out


def dersimonian_laird(effects: np.ndarray, ses: np.ndarray) -> dict:
    """Random-effects meta-analysis (DerSimonian–Laird)."""
    y = np.asarray(effects, dtype=float)
    se = np.asarray(ses, dtype=float)
    k = len(y)
    w = 1.0 / (se**2)
    ybar = np.sum(w * y) / np.sum(w)
    Q = float(np.sum(w * (y - ybar) ** 2))
    df = k - 1
    c = float(np.sum(w) - np.sum(w**2) / np.sum(w))
    tau2 = max(0.0, (Q - df) / c) if c > 0 else 0.0
    w_star = 1.0 / (se**2 + tau2)
    pooled = float(np.sum(w_star * y) / np.sum(w_star))
    se_p = float(1.0 / math.sqrt(np.sum(w_star)))
    lo = pooled - Z95 * se_p
    hi = pooled + Z95 * se_p
    z = pooled / se_p if se_p > 0 else 0.0
    p = float(2 * (1 - stats.norm.cdf(abs(z))))
    i2 = max(0.0, (Q - df) / Q) * 100.0 if Q > 0 else 0.0
    q_p = float(1 - stats.chi2.cdf(Q, df)) if df > 0 else 1.0
    weights = w_star / np.sum(w_star) * 100.0
    return {
        "pooled_effect": pooled,
        "pooled_SE": se_p,
        "lower_CI": lo,
        "upper_CI": hi,
        "p_value": p,
        "I2": i2,
        "tau2": tau2,
        "Q": Q,
        "Q_p": q_p,
        "n_studies": k,
        "weights_pct": weights,
    }


def direction_from_effect(e: float, se: float, alpha: float = 0.05) -> str:
    """Directional label with a weak-evidence 'no clear association' band."""
    z = e / se if se > 0 else 0.0
    p = 2 * (1 - stats.norm.cdf(abs(z)))
    if p >= alpha and abs(e) < 0.25:
        return "no_clear_association"
    if e > 0:
        return "increased_in_MetS"
    if e < 0:
        return "decreased_in_MetS"
    return "no_clear_association"


def arrow_symbol(direction: str) -> str:
    if direction == "increased_in_MetS":
        return "↑"
    if direction == "decreased_in_MetS":
        return "↓"
    return "—"


def load_real() -> tuple[pd.DataFrame, pd.DataFrame]:
    d1 = pd.read_csv(REAL_101, sep="\t")
    d2 = pd.read_csv(REAL_838, sep="\t")
    assert (d1["BioProject"] == "PRJNA1016449").all()
    assert (d2["BioProject"] == "PRJNA838778").all()
    return d1, d2


def load_real_three() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    d1, d2 = load_real()
    if not REAL_819.exists():
        raise FileNotFoundError(
            f"Real PRJNA819279 effect sizes not found: {REAL_819}\n"
            "Complete the single-study pipeline first."
        )
    d3 = pd.read_csv(REAL_819, sep="\t")
    assert (d3["BioProject"] == "PRJNA819279").all()
    return d1, d2, d3


def _majority_same_direction(effects: list[float]) -> bool:
    signs = [int(np.sign(e)) for e in effects if e != 0]
    if not signs:
        return False
    return max(signs.count(1), signs.count(-1)) >= 2


def select_taxa_three(
    d1: pd.DataFrame, d2: pd.DataFrame, d3: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Select genera present in all three real cohorts with valid estimates."""
    s1, s2, s3 = set(d1["taxon"]), set(d2["taxon"]), set(d3["taxon"])
    shared = sorted(s1 & s2 & s3)
    union = sorted(s1 | s2 | s3)

    harm_rows = []
    for t in union:
        in1, in2, in3 = t in s1, t in s2, t in s3
        if in1 and in2 and in3:
            match_type = "exact_shared_all_three"
        elif in1 and in2:
            match_type = "shared_1016449_838778"
        elif in1 and in3:
            match_type = "shared_1016449_819279"
        elif in2 and in3:
            match_type = "shared_838778_819279"
        else:
            match_type = "single_study_only"
        harm_rows.append(
            {
                "harmonized_taxon": t,
                "PRJNA1016449_taxon": t if in1 else "",
                "PRJNA838778_taxon": t if in2 else "",
                "PRJNA819279_taxon": t if in3 else "",
                "match_type": match_type,
                "selected_for_meta": "pending",
            }
        )
    harm = pd.DataFrame(harm_rows)

    m1, m2, m3 = d1.set_index("taxon"), d2.set_index("taxon"), d3.set_index("taxon")
    sel_rows = []
    for t in shared:
        a, b, c = m1.loc[t], m2.loc[t], m3.loc[t]
        prevs = [float(a["prevalence"]), float(b["prevalence"]), float(c["prevalence"])]
        effects = [float(a["effect_size"]), float(b["effect_size"]), float(c["effect_size"])]
        ses = [float(a["standard_error"]), float(b["standard_error"]), float(c["standard_error"])]
        prev_ok = all(p >= MIN_PREVALENCE for p in prevs)
        valid_est = all(
            np.isfinite(e) and np.isfinite(se) and se > 0 for e, se in zip(effects, ses)
        )
        same_majority = _majority_same_direction(effects)
        mean_abs = float(np.mean(np.abs(effects)))
        interpretable = (
            same_majority
            or mean_abs >= 0.2
            or min(prevs) >= 0.3
        )
        selected = bool(prev_ok and valid_est and interpretable)
        reason_parts = []
        if not prev_ok:
            reason_parts.append("low_prevalence")
        if not valid_est:
            reason_parts.append("invalid_estimate")
        if not interpretable:
            reason_parts.append("not_interpretable")
        if selected:
            if same_majority and mean_abs >= 0.3:
                reason_parts.append("majority_concordant_nontrivial")
            elif same_majority:
                reason_parts.append("majority_concordant_weak")
            elif mean_abs >= 0.2:
                reason_parts.append("heterogeneous_but_nontrivial")
            else:
                reason_parts.append("high_prevalence_shared")
        sel_rows.append(
            {
                "taxon": t,
                "taxonomic_level": "genus",
                "in_PRJNA1016449": True,
                "in_PRJNA838778": True,
                "in_PRJNA819279": True,
                "prevalence_1016449": prevs[0],
                "prevalence_838778": prevs[1],
                "prevalence_819279": prevs[2],
                "effect_1016449": effects[0],
                "effect_838778": effects[1],
                "effect_819279": effects[2],
                "SE_1016449": ses[0],
                "SE_838778": ses[1],
                "SE_819279": ses[2],
                "majority_same_direction": bool(same_majority),
                "mean_abs_effect": mean_abs,
                "selected": selected,
                "selection_reason": ";".join(reason_parts) if reason_parts else "selected",
            }
        )
    sel = pd.DataFrame(sel_rows).sort_values(
        ["selected", "majority_same_direction", "mean_abs_effect"],
        ascending=[False, False, False],
    )
    selected_set = set(sel.loc[sel["selected"], "taxon"])
    harm["selected_for_meta"] = harm["harmonized_taxon"].apply(
        lambda t: "yes" if t in selected_set else "no"
    )
    sel.to_csv(META / "taxon_selection.tsv", sep="\t", index=False)
    harm.to_csv(META / "taxon_harmonization.tsv", sep="\t", index=False)
    return sel, harm


def build_three_study(
    d1: pd.DataFrame,
    d2: pd.DataFrame,
    d3: pd.DataFrame,
    sel: pd.DataFrame,
    d3_status: str = "real",
) -> pd.DataFrame:
    taxa = list(sel.loc[sel["selected"], "taxon"])
    cols = [
        "BioProject",
        "taxon",
        "taxonomic_level",
        "MetS_n",
        "Control_n",
        "effect_size",
        "standard_error",
        "lower_CI",
        "upper_CI",
        "p_value",
        "FDR",
        "effect_direction",
        "prevalence",
        "analysis_method",
        "model",
        "covariates",
        "data_status",
    ]

    def subset(df: pd.DataFrame, status: str) -> pd.DataFrame:
        x = df[df["taxon"].isin(taxa)].copy()
        x["data_status"] = status
        # Harmonize Control_n column name already present
        rename = {"Control_n": "Control_n", "effect_direction": "effect_direction"}
        x = x.rename(columns=rename)
        keep = [c for c in cols if c in x.columns]
        return x[keep]

    a = subset(d1, "real")
    b = subset(d2, "real")
    c = subset(d3, d3_status)
    three = pd.concat([a, b, c], ignore_index=True)
    three = three.rename(columns={"Control_n": "control_n", "effect_direction": "direction"})
    three = three[
        [
            "BioProject",
            "taxon",
            "taxonomic_level",
            "MetS_n",
            "control_n",
            "effect_size",
            "standard_error",
            "lower_CI",
            "upper_CI",
            "p_value",
            "FDR",
            "direction",
            "prevalence",
            "analysis_method",
            "model",
            "covariates",
            "data_status",
        ]
    ]
    three.to_csv(META / "three_study_effects.tsv", sep="\t", index=False)
    return three


def run_meta(three: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    weight_rows = []
    for taxon, g in three.groupby("taxon"):
        g = g.set_index("BioProject")
        effects = g.loc[["PRJNA1016449", "PRJNA838778", "PRJNA819279"], "effect_size"].values
        ses = g.loc[["PRJNA1016449", "PRJNA838778", "PRJNA819279"], "standard_error"].values
        res = dersimonian_laird(effects, ses)
        direction = direction_from_effect(res["pooled_effect"], res["pooled_SE"])
        rows.append(
            {
                "taxon": taxon,
                "taxonomic_level": "genus",
                "pooled_effect": res["pooled_effect"],
                "pooled_SE": res["pooled_SE"],
                "lower_CI": res["lower_CI"],
                "upper_CI": res["upper_CI"],
                "p_value": res["p_value"],
                "I2": res["I2"],
                "tau2": res["tau2"],
                "Q": res["Q"],
                "Q_p": res["Q_p"],
                "n_studies": res["n_studies"],
                "direction": direction,
                "effect_PRJNA1016449": float(effects[0]),
                "effect_PRJNA838778": float(effects[1]),
                "effect_PRJNA819279": float(effects[2]),
                "weight_PRJNA1016449": float(res["weights_pct"][0]),
                "weight_PRJNA838778": float(res["weights_pct"][1]),
                "weight_PRJNA819279": float(res["weights_pct"][2]),
            }
        )
        for bp, w in zip(
            ["PRJNA1016449", "PRJNA838778", "PRJNA819279"], res["weights_pct"]
        ):
            weight_rows.append({"taxon": taxon, "BioProject": bp, "weight_pct": float(w)})

    meta = pd.DataFrame(rows)
    meta["FDR"] = bh_fdr(meta["p_value"].values)
    meta = meta.sort_values("p_value")
    # Column order
    meta = meta[
        [
            "taxon",
            "taxonomic_level",
            "pooled_effect",
            "pooled_SE",
            "lower_CI",
            "upper_CI",
            "p_value",
            "FDR",
            "I2",
            "tau2",
            "Q",
            "Q_p",
            "n_studies",
            "direction",
            "effect_PRJNA1016449",
            "effect_PRJNA838778",
            "effect_PRJNA819279",
            "weight_PRJNA1016449",
            "weight_PRJNA838778",
            "weight_PRJNA819279",
        ]
    ]
    meta.to_csv(META / "three_study_random_effects.tsv", sep="\t", index=False)

    het = meta[
        ["taxon", "pooled_effect", "lower_CI", "upper_CI", "I2", "tau2", "Q", "Q_p", "n_studies"]
    ].copy()
    het.to_csv(META / "heterogeneity_results.tsv", sep="\t", index=False)
    return meta, pd.DataFrame(weight_rows)


def classify_reproducibility(three: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, m in meta.iterrows():
        t = m["taxon"]
        g = three[three["taxon"] == t].set_index("BioProject")
        dirs = {}
        for bp in ["PRJNA1016449", "PRJNA838778", "PRJNA819279"]:
            e = float(g.loc[bp, "effect_size"])
            se = float(g.loc[bp, "standard_error"])
            dirs[bp] = direction_from_effect(e, se)
        signs = [np.sign(float(g.loc[bp, "effect_size"])) for bp in dirs]
        nonzero_signs = [s for s in signs if s != 0]
        same_all = len(set(nonzero_signs)) == 1 and len(nonzero_signs) == 3
        n_clear = sum(d != "no_clear_association" for d in dirs.values())
        pooled_clear = m["direction"] != "no_clear_association"
        i2 = float(m["I2"])
        abs_pool = abs(float(m["pooled_effect"]))
        fdr = float(m["FDR"])

        if same_all and pooled_clear and fdr < 0.2 and i2 < 60:
            cat = "strongly_concordant"
        elif same_all and pooled_clear:
            cat = "directionally_concordant"
        elif len(set(nonzero_signs)) > 1 or i2 >= 60:
            cat = "heterogeneous"
        else:
            cat = "no_clear_association"

        rows.append(
            {
                "taxon": t,
                "class": cat,
                "dir_PRJNA1016449": dirs["PRJNA1016449"],
                "dir_PRJNA838778": dirs["PRJNA838778"],
                "dir_PRJNA819279": dirs["PRJNA819279"],
                "dir_pooled": m["direction"],
                "symbol_PRJNA1016449": arrow_symbol(dirs["PRJNA1016449"]),
                "symbol_PRJNA838778": arrow_symbol(dirs["PRJNA838778"]),
                "symbol_PRJNA819279": arrow_symbol(dirs["PRJNA819279"]),
                "symbol_pooled": arrow_symbol(m["direction"]),
                "pooled_effect": m["pooled_effect"],
                "FDR": m["FDR"],
                "I2": m["I2"],
                "n_clear_cohorts": n_clear,
            }
        )
    clas = pd.DataFrame(rows).sort_values(["class", "FDR"])
    clas.to_csv(META / "reproducibility_classification.tsv", sep="\t", index=False)
    return clas


def rank_top_taxa(meta: pd.DataFrame, clas: pd.DataFrame) -> pd.DataFrame:
    c = clas.set_index("taxon")
    rows = []
    class_score = {
        "strongly_concordant": 4,
        "directionally_concordant": 3,
        "heterogeneous": 1,
        "no_clear_association": 0,
    }
    for _, m in meta.iterrows():
        t = m["taxon"]
        cs = class_score[c.loc[t, "class"]]
        # Evidence score: consistency, magnitude, precision, FDR, low heterogeneity
        mag = abs(float(m["pooled_effect"]))
        precision = 1.0 / (1.0 + float(m["pooled_SE"]))
        fdr_score = 1.0 - float(m["FDR"])
        het_score = 1.0 - min(float(m["I2"]), 100) / 100.0
        n_sup = int(c.loc[t, "n_clear_cohorts"])
        score = (
            2.5 * cs
            + 1.5 * mag
            + 1.0 * precision
            + 1.2 * fdr_score
            + 1.0 * het_score
            + 0.5 * n_sup
        )
        rows.append(
            {
                "taxon": t,
                "rank_score": score,
                "class": c.loc[t, "class"],
                "pooled_effect": m["pooled_effect"],
                "lower_CI": m["lower_CI"],
                "upper_CI": m["upper_CI"],
                "p_value": m["p_value"],
                "FDR": m["FDR"],
                "I2": m["I2"],
                "n_supporting_cohorts": n_sup,
                "consistency_score": cs,
                "direction": m["direction"],
            }
        )
    top = pd.DataFrame(rows).sort_values("rank_score", ascending=False)
    top.insert(0, "rank", range(1, len(top) + 1))
    top.to_csv(META / "top_reproducible_taxa.tsv", sep="\t", index=False)
    return top


def main() -> None:
    META.mkdir(parents=True, exist_ok=True)
    d1, d2, d3 = load_real_three()
    sel, _harm = select_taxa_three(d1, d2, d3)
    print(f"{int(sel['selected'].sum())} of {len(sel)} shared genera entered the meta-analysis")
    three = build_three_study(d1, d2, d3, sel)
    meta, _weights = run_meta(three)
    clas = classify_reproducibility(three, meta)
    rank_top_taxa(meta, clas)
    print(f"Tables written to {META}")


if __name__ == "__main__":
    main()
