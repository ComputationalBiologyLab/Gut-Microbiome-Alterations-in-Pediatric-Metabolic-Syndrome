#!/usr/bin/env python3
"""Random-effects meta-analysis at phylum level from ANCOM-BC2 study-level effects."""

from __future__ import annotations

import math
from pathlib import Path

import os
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(os.environ.get("MS_ROOT", ".")).resolve()  # directory with the <BioProject>_analysis folders
META = ROOT / "meta_analysis"
Z95 = 1.959963984540054


def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    adj = np.empty(n, dtype=float)
    prev = 1.0
    for i in range(n - 1, -1, -1):
        val = ranked[i] * n / (i + 1)
        prev = min(prev, val)
        adj[i] = prev
    out = np.empty(n, dtype=float)
    out[order] = np.clip(adj, 0, 1)
    return out


def dersimonian_laird(effects: np.ndarray, ses: np.ndarray) -> dict:
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
    }


def main() -> None:
    src = META / "phylum_study_level_effects_all.tsv"
    if not src.exists():
        raise SystemExit(f"Missing {src}; run run_phylum_ancombc2_all_studies.R first")

    t = pd.read_csv(src, sep="\t")
    # Shared phyla present in all three with valid estimates
    shared = sorted(
        set(t.groupby("taxon").filter(lambda g: len(g) == 3)["taxon"].unique())
    )
    rows = []
    for phylum in shared:
        g = t[t.taxon == phylum].set_index("BioProject")
        effects = g.loc[["PRJNA1016449", "PRJNA838778", "PRJNA819279"], "effect_size"].values.astype(float)
        ses = g.loc[["PRJNA1016449", "PRJNA838778", "PRJNA819279"], "standard_error"].values.astype(float)
        if not (np.all(np.isfinite(effects)) and np.all(ses > 0)):
            continue
        res = dersimonian_laird(effects, ses)
        rows.append(
            {
                "taxon": phylum,
                "taxonomic_level": "phylum",
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
                "effect_PRJNA1016449": effects[0],
                "effect_PRJNA838778": effects[1],
                "effect_PRJNA819279": effects[2],
            }
        )
    meta = pd.DataFrame(rows).sort_values("p_value")
    meta["FDR"] = bh_fdr(meta["p_value"].values)
    meta.to_csv(META / "phylum_three_study_random_effects.tsv", sep="\t", index=False)

    print("=== PHYLUM-LEVEL META-ANALYSIS ===")
    print(f"Phyla tested: {len(meta)}")
    print(meta[["taxon", "pooled_effect", "p_value", "FDR", "I2", "lower_CI", "upper_CI"]].to_string(index=False))
    print(f"\nFDR < 0.10: {(meta.FDR < 0.1).sum()}")
    print(f"FDR < 0.20: {(meta.FDR < 0.2).sum()}")
    print(f"Nominal p < 0.05: {(meta.p_value < 0.05).sum()}")

    # Single-study phylum FDR
    print("\n=== SINGLE-STUDY PHYLUM FDR < 0.1 ===")
    for bp in ["PRJNA1016449", "PRJNA838778", "PRJNA819279"]:
        sub = t[t.BioProject == bp]
        sig = sub[sub.FDR < 0.1]
        print(f"{bp}: {list(sig.taxon)} (n={len(sig)})")


if __name__ == "__main__":
    main()
