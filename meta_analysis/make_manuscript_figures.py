#!/usr/bin/env python3
"""Build manuscript figures: single-study collages (native re-render) + composites."""
import os, numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# MS_ROOT: directory that contains the <BioProject>_analysis folders and meta_analysis/ (default: current directory)
BASE = os.path.abspath(os.environ.get("MS_ROOT", "."))
OUT = os.path.join(BASE, "figures", "manuscript")
os.makedirs(OUT, exist_ok=True)

COHORTS = [
    ("PRJNA1016449", "USA", "V3–V4"),
    ("PRJNA838778",  "China", "V4"),
    ("PRJNA819279",  "Mexico", "V3–V4"),
]
C_CTRL, C_METS = "#2C7FB8", "#D95F02"
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8,
    "axes.labelsize": 8.5, "axes.titlesize": 9, "xtick.labelsize": 8,
    "ytick.labelsize": 8, "legend.fontsize": 8, "axes.linewidth": 0.8,
    "xtick.major.width": 0.8, "ytick.major.width": 0.8, "savefig.dpi": 600,
})

def alpha(bp):
    return pd.read_csv(f"{BASE}/{bp}_analysis/results/alpha_diversity_rarefied.tsv", sep="\t")

def alpha_stats(bp):
    t = pd.read_csv(f"{BASE}/{bp}_analysis/12_tables/Table3_alpha_diversity.tsv", sep="\t")
    return t.set_index("metric")

def permanova(bp):
    t = pd.read_csv(f"{BASE}/{bp}_analysis/12_tables/Table4_beta_diversity_permanova.tsv", sep="\t")
    return t.set_index("metric").loc["braycurtis"]

def pcoa(bp):
    """Classical MDS from the Bray-Curtis distance matrix."""
    d = pd.read_csv(f"{BASE}/{bp}_analysis/07_diversity/bray_curtis_export/distance-matrix.tsv",
                    sep="\t", index_col=0)
    D = d.values.astype(float); n = D.shape[0]
    J = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * J @ (D ** 2) @ J
    w, V = np.linalg.eigh(B)
    idx = np.argsort(w)[::-1]; w, V = w[idx], V[:, idx]
    pos = w > 0
    coords = V[:, pos] * np.sqrt(w[pos])
    expl = w[pos] / w[pos].sum() * 100
    return pd.DataFrame(coords[:, :2], index=d.index, columns=["PC1", "PC2"]), expl[:2]

def effects(bp):
    return pd.read_csv(f"{BASE}/{bp}_analysis/results/study_level_effect_sizes.tsv", sep="\t")

def panel_tag(ax, letter):
    ax.annotate(letter, xy=(0, 1), xycoords="axes fraction", xytext=(-30, 5),
                textcoords="offset points", fontsize=11, fontweight="bold",
                va="bottom", ha="left", annotation_clip=False)

# ---------------------------------------------------------------- Figure 2
def figure2():
    fig, axes = plt.subplots(2, 3, figsize=(7.4, 5.4), constrained_layout=True)
    letters = iter("ABCDEF")
    for j, (bp, country, region) in enumerate(COHORTS):
        a = alpha(bp); st = alpha_stats(bp)
        ctrl = a.loc[a.mets_status == "Control", "shannon"].values
        mets = a.loc[a.mets_status == "MetS", "shannon"].values
        ax = axes[0, j]
        bpx = ax.boxplot([ctrl, mets], positions=[0, 1], widths=0.55, patch_artist=True,
                         showfliers=False, medianprops=dict(color="black", lw=1.2),
                         whiskerprops=dict(lw=0.8), capprops=dict(lw=0.8),
                         boxprops=dict(lw=0.8))
        for patch, c in zip(bpx["boxes"], [C_CTRL, C_METS]):
            patch.set_facecolor(c); patch.set_alpha(0.30)
        rng = np.random.default_rng(0)
        for k, (v, c) in enumerate([(ctrl, C_CTRL), (mets, C_METS)]):
            ax.scatter(rng.normal(k, 0.075, len(v)), v, s=9, color=c,
                       edgecolor="white", linewidth=0.3, zorder=3, alpha=0.9)
        ax.set_xticks([0, 1])
        ax.set_xticklabels([f"Control\n(n={len(ctrl)})", f"MetS\n(n={len(mets)})"])
        ax.set_ylabel("Shannon index" if j == 0 else "")
        lo, hi = ax.get_ylim(); ax.set_ylim(lo, hi + 0.22 * (hi - lo))
        ax.text(0.5, 0.985, f"FDR = {float(st.loc['shannon','FDR']):.3f}",
                transform=ax.transAxes, ha="center", va="top", fontsize=7.5,
                bbox=dict(fc="white", ec="0.8", lw=0.5, pad=1.6))
        ax.set_title(f"{bp}\n{country} | 16S {region}", fontsize=8.5, fontweight="bold", pad=5)
        panel_tag(ax, next(letters))
        ax.spines[["top", "right"]].set_visible(False)

    for j, (bp, country, region) in enumerate(COHORTS):
        co, expl = pcoa(bp); a = alpha(bp).set_index("sample_id")
        grp = a.reindex(co.index)["mets_status"]; pm = permanova(bp)
        ax = axes[1, j]
        for g, c in [("Control", C_CTRL), ("MetS", C_METS)]:
            m = (grp == g).values
            ax.scatter(co.loc[m, "PC1"], co.loc[m, "PC2"], s=16, color=c,
                       edgecolor="white", linewidth=0.4, alpha=0.9, label=g)
        ax.set_xlabel(f"PCo1 ({expl[0]:.1f}%)"); ax.set_ylabel(f"PCo2 ({expl[1]:.1f}%)")
        lo, hi = ax.get_ylim(); ax.set_ylim(lo, hi + 0.26 * (hi - lo))
        ax.text(0.03, 0.985, f"PERMANOVA R² = {float(pm['R2']):.3f}\np = {float(pm['p_value']):.3f}",
                transform=ax.transAxes, ha="left", va="top", fontsize=7.0,
                bbox=dict(fc="white", ec="0.8", lw=0.5, pad=1.6))
        panel_tag(ax, next(letters))
        ax.spines[["top", "right"]].set_visible(False)
    h = [Line2D([], [], marker="o", ls="", color=C_CTRL, label="Control (non-MetS)"),
         Line2D([], [], marker="o", ls="", color=C_METS, label="MetS")]
    fig.legend(handles=h, ncol=2, frameon=False, loc="outside lower center")
    for ext in ("png", "pdf"):
        fig.savefig(f"{OUT}/Figure2_single_study_community.{ext}", bbox_inches="tight")
    plt.close(fig); print("Figure 2 written")

# ---------------------------------------------------------------- Figure 3
def figure3():
    fig, axes = plt.subplots(2, 3, figsize=(7.4, 6.4), constrained_layout=True)
    letters = iter("ABCDEF")
    for j, (bp, country, region) in enumerate(COHORTS):
        e = effects(bp).copy()
        e["nlf"] = -np.log10(e["FDR"].clip(lower=1e-12))
        ax = axes[0, j]
        sig = e["FDR"] < 0.1
        ax.scatter(e.loc[~sig, "effect_size"], e.loc[~sig, "nlf"], s=10, color="0.72",
                   edgecolor="none", alpha=0.85)
        up = sig & (e.effect_size > 0); dn = sig & (e.effect_size < 0)
        ax.scatter(e.loc[up, "effect_size"], e.loc[up, "nlf"], s=18, color=C_METS,
                   edgecolor="white", linewidth=0.35, zorder=3)
        ax.scatter(e.loc[dn, "effect_size"], e.loc[dn, "nlf"], s=18, color=C_CTRL,
                   edgecolor="white", linewidth=0.35, zorder=3)
        ax.axhline(-np.log10(0.1), ls="--", lw=0.7, color="0.35")
        ax.axvline(0, ls=":", lw=0.7, color="0.55")
        ax.set_ylim(-0.05, max(e["nlf"].max(), 1.05) * 1.30)
        ax.set_xlabel("ANCOM-BC2 log fold-change")
        ax.set_ylabel("−log$_{10}$ FDR" if j == 0 else "")
        ax.set_title(f"{bp}\n{country}", fontsize=8.5, fontweight="bold", pad=5)
        ax.text(0.03, 0.985, f"{int(sig.sum())}/{len(e)} genera at FDR < 0.1",
                transform=ax.transAxes, ha="left", va="top", fontsize=7.0,
                bbox=dict(fc="white", ec="0.8", lw=0.5, pad=1.6))
        ax.text(0.97, -np.log10(0.1), "FDR 0.1", transform=ax.get_yaxis_transform(),
                ha="right", va="bottom", fontsize=6.6, color="0.35")
        ax.locator_params(axis="x", nbins=3)
        ax.tick_params(axis="x", labelsize=7)
        panel_tag(ax, next(letters))
        ax.spines[["top", "right"]].set_visible(False)

    for j, (bp, country, region) in enumerate(COHORTS):
        e = effects(bp).copy()
        top = e.reindex(e["effect_size"].abs().sort_values(ascending=False).index).head(10)
        top = top.sort_values("effect_size")
        ax = axes[1, j]
        cols = [C_METS if v > 0 else C_CTRL for v in top.effect_size]
        y = np.arange(len(top))
        ax.barh(y, top.effect_size, color=cols, alpha=0.88, height=0.68,
                edgecolor="white", linewidth=0.4)
        ax.errorbar(top.effect_size, y,
                    xerr=[top.effect_size - top.lower_CI, top.upper_CI - top.effect_size],
                    fmt="none", ecolor="0.25", elinewidth=0.7, capsize=1.6)
        ax.set_yticks(y)
        ax.set_yticklabels([t if len(t) <= 27 else t[:25] + "…" for t in top.taxon], fontsize=6.4)
        lo, hi = float(top.lower_CI.min()), float(top.upper_CI.max())
        span = hi - lo
        ax.set_xlim(lo - 0.16 * span, hi + 0.16 * span)
        for k, (f, ev, u, l) in enumerate(zip(top.FDR, top.effect_size, top.upper_CI, top.lower_CI)):
            xx = (u + 0.03 * span) if ev > 0 else (l - 0.03 * span)
            if f < 0.1:
                ax.text(xx, k, "*", ha="left" if ev > 0 else "right", va="center",
                        fontsize=10, color="0.15")
        ax.axvline(0, color="0.3", lw=0.8)
        ax.set_xlabel("ANCOM-BC2 log fold-change")
        ax.locator_params(axis="x", nbins=3)
        ax.tick_params(axis="x", labelsize=7)
        panel_tag(ax, next(letters))
        ax.spines[["top", "right"]].set_visible(False)
    h = [Line2D([], [], marker="s", ls="", color=C_METS, label="Increased in MetS"),
         Line2D([], [], marker="s", ls="", color=C_CTRL, label="Decreased in MetS"),
         Line2D([], [], marker="$*$", ls="", color="0.15", label="FDR < 0.1")]
    fig.legend(handles=h, ncol=3, frameon=False, loc="outside lower center")
    for ext in ("png", "pdf"):
        fig.savefig(f"{OUT}/Figure3_single_study_differential_abundance.{ext}", bbox_inches="tight")
    plt.close(fig); print("Figure 3 written")

figure2(); figure3()


# ---------------------------------------------------------------- Figure 4
META = os.path.join(BASE, "meta_analysis")
ORDER = ["[Eubacterium]_xylanophilum_group", "Alistipes", "Terrisporobacter", "Dialister",
         "UCG-002", "Odoribacter", "CAG-56", "[Eubacterium]_ruminantium_group",
         "Faecalibacterium", "[Eubacterium]_eligens_group", "Romboutsia"]
BPS = ["PRJNA1016449", "PRJNA838778", "PRJNA819279"]
LBL = {"PRJNA1016449": "PRJNA1016449 (USA)", "PRJNA838778": "PRJNA838778 (China)",
       "PRJNA819279": "PRJNA819279 (Mexico)"}

def meta_tbl():
    m = pd.read_csv(f"{META}/three_study_random_effects.tsv", sep="\t").set_index("taxon")
    return m

def study_se(bp, taxon):
    e = effects(bp).set_index("taxon")
    return float(e.loc[taxon, "standard_error"])

RANK = {"f__": "family", "o__": "order", "c__": "class", "p__": "phylum", "g__": "genus"}

def pretty(t):
    for pre, rank in RANK.items():
        if t.startswith(pre):
            return f"{t[len(pre):]} ({rank})"
    return t.replace("_group", " group").replace("_", " ")

def figure4():
    m = meta_tbl()
    se_cache = {bp: effects(bp).set_index("taxon")["standard_error"] for bp in BPS}
    rows = []
    for t in ORDER:
        for bp in BPS:
            rows.append((t, LBL[bp], float(m.loc[t, f"effect_{bp}"]), float(se_cache[bp][t]), False))
        rows.append((t, "Pooled (random effects)", float(m.loc[t, "pooled_effect"]),
                     float(m.loc[t, "pooled_SE"]), True))
    n = len(rows)
    fig = plt.figure(figsize=(7.4, 9.2))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 0.86], wspace=0.02)
    ax = fig.add_subplot(gs[0, 0]); axt = fig.add_subplot(gs[0, 1]); axt.axis("off")
    ypos = {}
    y = 0.0
    for i, (t, lab, e, se, pooled) in enumerate(rows):
        if i and rows[i - 1][0] != t:
            y += 1.15
        ypos[i] = y; y += 1.0
    ymax = y
    for i, (t, lab, e, se, pooled) in enumerate(rows):
        yy = ymax - ypos[i]
        lo, hi = e - 1.96 * se, e + 1.96 * se
        col = "#1F4E79" if pooled else "0.35"
        ax.plot([lo, hi], [yy, yy], color=col, lw=1.5 if pooled else 1.0,
                solid_capstyle="butt", zorder=2)
        if pooled:
            ax.plot([e], [yy], marker="D", ms=6, color=col, zorder=3)
        else:
            ax.plot([e], [yy], marker="s", ms=4.2,
                    color=C_METS if e > 0 else C_CTRL, zorder=3)
        ax.text(-0.012, yy, lab, transform=ax.get_yaxis_transform(), ha="right",
                va="center", fontsize=6.6,
                fontweight="bold" if pooled else "normal",
                color="#1F4E79" if pooled else "0.25")
        axt.text(0.02, yy, f"{e:+.3f}", transform=axt.get_yaxis_transform(),
                 ha="left", va="center", fontsize=6.6,
                 fontweight="bold" if pooled else "normal")
        axt.text(0.30, yy, f"({lo:+.2f} to {hi:+.2f})", transform=axt.get_yaxis_transform(),
                 ha="left", va="center", fontsize=6.6,
                 fontweight="bold" if pooled else "normal")
        if pooled:
            axt.text(0.74, yy, f"{float(m.loc[t,'I2']):.1f}", transform=axt.get_yaxis_transform(),
                     ha="right", va="center", fontsize=6.6, fontweight="bold")
            axt.text(0.99, yy, f"{float(m.loc[t,'FDR']):.3f}", transform=axt.get_yaxis_transform(),
                     ha="right", va="center", fontsize=6.6, fontweight="bold")
    for t in ORDER:
        idx = [i for i, r in enumerate(rows) if r[0] == t]
        yy = ymax - ypos[idx[0]] + 0.95
        ax.text(-0.012, yy, pretty(t), transform=ax.get_yaxis_transform(), ha="right",
                va="center", fontsize=7.4, fontweight="bold", style="italic", color="black")
    ax.axvline(0, color="0.3", lw=0.9, ls="--")
    ax.set_ylim(-0.8, ymax + 1.8); ax.set_yticks([])
    ax.set_xlim(-2.1, 2.1)
    ax.set_xlabel("ANCOM-BC2 log fold-change (MetS vs Control), 95% CI")
    ax.spines[["top", "right", "left"]].set_visible(False)
    axt.set_ylim(*ax.get_ylim())
    hdr = ymax + 1.15
    axt.text(0.02, hdr, "Effect", transform=axt.get_yaxis_transform(), fontsize=7,
             fontweight="bold", ha="left", va="center")
    axt.text(0.30, hdr, "95% CI", transform=axt.get_yaxis_transform(), fontsize=7,
             fontweight="bold", ha="left", va="center")
    axt.text(0.74, hdr, "I² (%)", transform=axt.get_yaxis_transform(), fontsize=7,
             fontweight="bold", ha="right", va="center")
    axt.text(0.99, hdr, "FDR", transform=axt.get_yaxis_transform(), fontsize=7,
             fontweight="bold", ha="right", va="center")
    ax.text(0.5, hdr, "Decreased in MetS  ←     → Increased in MetS",
            transform=ax.get_yaxis_transform(), fontsize=7, ha="center", va="center",
            color="0.3")
    for ext in ("png", "pdf"):
        fig.savefig(f"{OUT}/Figure4_forest_concordant_genera.{ext}", bbox_inches="tight")
    plt.close(fig); print("Figure 4 (forest) written")

figure4()


# ---------------------------------------------------------------- Figure 5
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import Rectangle

CLASS_ORDER = ["directionally_concordant", "heterogeneous", "no_clear_association"]
CLASS_LAB = {"directionally_concordant": "Directionally\nconcordant (11)",
             "heterogeneous": "Heterogeneous (45)",
             "no_clear_association": "No clear\nassociation (5)"}
CLASS_COL = {"directionally_concordant": "#1F4E79", "heterogeneous": "#B08A2E",
             "no_clear_association": "0.55"}

def figure5():
    m = meta_tbl()
    cl = pd.read_csv(f"{META}/reproducibility_classification.tsv", sep="\t").set_index("taxon")
    m = m.join(cl["class"])
    m["absE"] = m.pooled_effect.abs()
    fig = plt.figure(figsize=(7.4, 9.8))
    gs = fig.add_gridspec(2, 3, width_ratios=[1.42, 0.74, 0.74],
                          height_ratios=[1.42, 1.0], wspace=0.60, hspace=0.26)
    axA = fig.add_subplot(gs[:, 0])

    blocks, ylabels = [], []
    for c in CLASS_ORDER:
        sub = m[m["class"] == c].sort_values("pooled_effect")
        blocks.append((c, sub)); ylabels += list(sub.index)
    mat = np.vstack([np.array([[float(r[f"effect_{bp}"]) for bp in BPS]
                               for _, r in sub.iterrows()]) for _, sub in blocks])
    vmax = np.nanpercentile(np.abs(mat), 97)
    im = axA.imshow(mat, aspect="auto", cmap="RdBu_r",
                    norm=TwoSlopeNorm(vcenter=0, vmin=-vmax, vmax=vmax))
    axA.set_xticks(range(3))
    axA.set_xticklabels(["USA", "China", "Mexico"], fontsize=7.5)
    axA.set_yticks(range(len(ylabels)))
    axA.set_yticklabels([pretty(t) for t in ylabels], fontsize=4.9, style="italic")
    axA.tick_params(axis="y", length=1.5, pad=1.5)
    for lab, t in zip(axA.get_yticklabels(), ylabels):
        if t in ORDER:
            lab.set_fontweight("bold"); lab.set_color("#1F4E79")
    start = 0
    for c, sub in blocks:
        end = start + len(sub)
        axA.add_patch(Rectangle((-0.5, start - 0.5), 3, len(sub), fill=False,
                                ec=CLASS_COL[c], lw=1.6, zorder=5))
        axA.annotate(CLASS_LAB[c], xy=(2.62, (start + end - 1) / 2), rotation=90,
                     va="center", ha="center", fontsize=6.6, fontweight="bold",
                     color=CLASS_COL[c], annotation_clip=False)
        start = end
    cb = fig.colorbar(im, ax=axA, orientation="horizontal", pad=0.045, fraction=0.032)
    cb.set_label("Cohort-level ANCOM-BC2 log fold-change\n"
                 "16/61 genera (26.2%) share direction in all three cohorts;\n"
                 "25% expected by chance (exact binomial p = 0.88)", fontsize=6.8)
    cb.ax.tick_params(labelsize=6.5)
    axA.annotate("A", xy=(0, 1), xycoords="axes fraction", xytext=(-62, 8),
                 textcoords="offset points", fontsize=11, fontweight="bold",
                 annotation_clip=False)

    axB = fig.add_subplot(gs[0, 1])
    for c in CLASS_ORDER:
        s2 = m[m["class"] == c]
        axB.scatter(s2.pooled_effect, -np.log10(s2.FDR), s=16, color=CLASS_COL[c],
                    edgecolor="white", linewidth=0.3, alpha=0.9)
    axB.axhline(-np.log10(0.05), ls="--", lw=0.8, color=C_METS)
    axB.axhline(-np.log10(0.20), ls=":", lw=0.8, color="0.4")
    axB.text(0.98, -np.log10(0.05), "FDR 0.05 ", transform=axB.get_yaxis_transform(),
             ha="right", va="bottom", fontsize=6.4, color=C_METS, zorder=6)
    axB.text(0.98, -np.log10(0.20), " FDR 0.20", transform=axB.get_yaxis_transform(),
             ha="right", va="bottom", fontsize=6.4, color="0.4")
    axB.axvline(0, ls=":", lw=0.7, color="0.6")
    axB.text(0.03, 0.43, "3/61 at nominal p < 0.05\n(3.05 expected under null);\n0 at FDR < 0.20",
             transform=axB.transAxes, ha="left", va="top", fontsize=6.4,
             bbox=dict(fc="white", ec="0.8", lw=0.5, pad=1.6))
    axB.set_ylim(0, 1.45); axB.set_xlabel("Pooled log fold-change")
    axB.set_ylabel("−log$_{10}$ FDR")
    axB.spines[["top", "right"]].set_visible(False)
    axB.annotate("B", xy=(0, 1), xycoords="axes fraction", xytext=(-34, 6),
                 textcoords="offset points", fontsize=11, fontweight="bold", annotation_clip=False)

    axC = fig.add_subplot(gs[0, 2])
    lim = max(m.absE.max(), m.pooled_SE.max()) * 1.08
    axC.plot([0, lim], [0, lim], ls="--", lw=0.8, color="0.45")
    for c in CLASS_ORDER:
        s2 = m[m["class"] == c]
        axC.scatter(s2.absE, s2.pooled_SE, s=16, color=CLASS_COL[c], edgecolor="white",
                    linewidth=0.3, alpha=0.9)
    axC.set_xlim(0, lim); axC.set_ylim(0, lim)
    axC.set_xlabel("|Pooled effect|"); axC.set_ylabel("Pooled standard error")
    n_above = int((m.pooled_SE > m.absE).sum())
    axC.text(0.96, 0.06, f"SE > |effect|\nin {n_above}/61 genera", transform=axC.transAxes,
             ha="right", va="bottom", fontsize=6.8,
             bbox=dict(fc="white", ec="0.8", lw=0.5, pad=1.6))
    axC.spines[["top", "right"]].set_visible(False)
    axC.annotate("C", xy=(0, 1), xycoords="axes fraction", xytext=(-34, 6),
                 textcoords="offset points", fontsize=11, fontweight="bold", annotation_clip=False)

    axD = fig.add_subplot(gs[1, 1])
    axD.hist(m.I2, bins=np.arange(0, 105, 7.5), color="0.72", edgecolor="white", linewidth=0.6)
    axD.axvline(float(np.median(m.I2)), color=C_METS, lw=1.3)
    axD.axvline(60, color="0.3", ls="--", lw=0.9)
    axD.text(float(np.median(m.I2)) + 2.5, axD.get_ylim()[1] * 0.92,
             f"median {np.median(m.I2):.1f}%", fontsize=6.6, color=C_METS)
    axD.text(61.5, axD.get_ylim()[1] * 0.58, "I² = 60%\n(heterogeneity\nrule)", fontsize=6.2, color="0.3")
    axD.set_xlabel("I² (%)"); axD.set_ylabel("Number of genera")
    axD.spines[["top", "right"]].set_visible(False)
    axD.annotate("D", xy=(0, 1), xycoords="axes fraction", xytext=(-34, 6),
                 textcoords="offset points", fontsize=11, fontweight="bold", annotation_clip=False)

    axE = fig.add_subplot(gs[1, 2])
    pairs = [("PRJNA1016449", "PRJNA838778", "o", "USA–China"),
             ("PRJNA1016449", "PRJNA819279", "^", "USA–Mexico"),
             ("PRJNA838778", "PRJNA819279", "s", "China–Mexico")]
    agree = 0; tot = 0
    for a, b, mk, lab in pairs:
        x = m[f"effect_{a}"].astype(float); y = m[f"effect_{b}"].astype(float)
        agree += int((np.sign(x) == np.sign(y)).sum()); tot += len(x)
        axE.scatter(x, y, s=13, marker=mk, alpha=0.65, linewidth=0.3,
                    edgecolor="white", label=lab)
    axE.axhline(0, lw=0.7, color="0.55"); axE.axvline(0, lw=0.7, color="0.55")
    axE.set_xlabel("Log fold-change, cohort 1"); axE.set_ylabel("Log fold-change, cohort 2")
    axE.text(0.03, 0.97, f"Same direction in {agree}/{tot}\npairings ({100*agree/tot:.1f}%);\n50% expected by chance",
             transform=axE.transAxes, ha="left", va="top", fontsize=6.6,
             bbox=dict(fc="white", ec="0.8", lw=0.5, pad=1.6))
    axE.legend(frameon=False, fontsize=6.0, loc="lower right", handletextpad=0.1,
               borderpad=0.1, labelspacing=0.25)
    axE.spines[["top", "right"]].set_visible(False)
    axE.annotate("E", xy=(0, 1), xycoords="axes fraction", xytext=(-34, 6),
                 textcoords="offset points", fontsize=11, fontweight="bold", annotation_clip=False)

    h = [Line2D([], [], marker="o", ls="", color=CLASS_COL[c],
                label=CLASS_LAB[c].replace("\n", " ")) for c in CLASS_ORDER]
    fig.legend(handles=h, ncol=3, frameon=False, loc="lower center",
               bbox_to_anchor=(0.5, -0.005), fontsize=7.2)
    for ext in ("png", "pdf"):
        fig.savefig(f"{OUT}/Figure5_reproducibility_heterogeneity.{ext}", bbox_inches="tight")
    plt.close(fig); print("Figure 5 (reproducibility) written")

figure5()


# ---------------------------------------------------------------- Figure 1
from PIL import Image
import matplotlib.image as mpimg

def figure1():
    fp = os.path.join(BASE, "figures")
    prisma = Image.open(f"{fp}/PRISMA_flow_diagram.png").convert("RGB")
    flow   = Image.open(f"{fp}/08_meta_analysis_workflow.png").convert("RGB")
    pop    = Image.open(f"{fp}/09_study_population.png").convert("RGB")
    fig = plt.figure(figsize=(7.4, 4.6))
    gs = fig.add_gridspec(2, 2, width_ratios=[0.98, 1.12],
                          height_ratios=[1.0, 0.86], wspace=0.03, hspace=0.10)
    axP = fig.add_subplot(gs[:, 0]); axW = fig.add_subplot(gs[0, 1]); axS = fig.add_subplot(gs[1, 1])
    for ax, im, tag in [(axP, prisma, "A"), (axW, flow, "B"), (axS, pop, "C")]:
        ax.imshow(np.asarray(im)); ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_visible(False)
        ax.annotate(tag, xy=(0, 1), xycoords="axes fraction", xytext=(2, 4),
                    textcoords="offset points", fontsize=12, fontweight="bold",
                    va="bottom", ha="left", annotation_clip=False)
    for ext in ("png", "pdf"):
        fig.savefig(f"{OUT}/Figure1_selection_and_design.{ext}", bbox_inches="tight", dpi=600)
    plt.close(fig); print("Figure 1 written")

figure1()


# ---------------------------------------------------------------- Figure 4b: non-replication
def figure_nonreplication():
    E = {bp: effects(bp).set_index("taxon") for bp in BPS}
    sig = {bp: set(E[bp].index[E[bp].FDR < 0.1]) for bp in BPS}
    allsig = sorted(set().union(*sig.values()))
    order = []
    for bp in BPS:
        order += sorted([t for t in allsig if t in sig[bp] and t not in order])
    mat = np.full((len(order), 3), np.nan)
    star = np.zeros_like(mat, dtype=bool)
    for i, t in enumerate(order):
        for j, bp in enumerate(BPS):
            if t in E[bp].index:
                mat[i, j] = float(E[bp].loc[t, "effect_size"])
                star[i, j] = t in sig[bp]

    fig = plt.figure(figsize=(7.4, 6.4))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.18, 0.86], height_ratios=[1.0, 0.82],
                          wspace=0.55, hspace=0.42)
    axA = fig.add_subplot(gs[:, 0])
    vmax = np.nanpercentile(np.abs(mat), 96)
    masked = np.ma.masked_invalid(mat)
    cmap = plt.get_cmap("RdBu_r").copy(); cmap.set_bad("0.90")
    im = axA.imshow(masked, aspect="auto", cmap=cmap,
                    norm=TwoSlopeNorm(vcenter=0, vmin=-vmax, vmax=vmax))
    for i in range(len(order)):
        for j in range(3):
            if star[i, j]:
                axA.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False,
                                        ec="black", lw=1.5, zorder=5))
            if np.isnan(mat[i, j]):
                axA.text(j, i, "·", ha="center", va="center", fontsize=7, color="0.45")
    axA.set_xticks(range(3)); axA.set_xticklabels(["USA", "China", "Mexico"], fontsize=7.5)
    axA.set_yticks(range(len(order)))
    axA.set_yticklabels([pretty(t) for t in order], fontsize=6.0, style="italic")
    axA.tick_params(axis="y", length=1.5, pad=1.5)
    start = 0
    for bp, lab in zip(BPS, ["USA", "China", "Mexico"]):
        n = len([t for t in order if t in sig[bp] and order.index(t) >= start])
        n = sum(1 for t in order[start:start + 40] if t in sig[bp])
        n = len([t for t in order if t in sig[bp]])
        axA.add_patch(Rectangle((-0.5, start - 0.5), 3, n, fill=False, ec="0.25",
                                lw=1.0, ls=":", zorder=4))
        axA.annotate(f"discovered\nin {lab} ({n})", xy=(2.60, start + n / 2 - 0.5), rotation=90,
                     va="center", ha="center", fontsize=6.4, fontweight="bold",
                     color="0.25", annotation_clip=False)
        start += n
    cb = fig.colorbar(im, ax=axA, orientation="horizontal", pad=0.055, fraction=0.035)
    cb.set_label("ANCOM-BC2 log fold-change\nBlack outline: FDR < 0.1 in that cohort. "
                 "Grey: genus not detected.", fontsize=6.8)
    cb.ax.tick_params(labelsize=6.5)
    axA.annotate("A", xy=(0, 1), xycoords="axes fraction", xytext=(-70, 8),
                 textcoords="offset points", fontsize=11, fontweight="bold",
                 annotation_clip=False)

    keep = [t for t in order if all(t in E[bp].index for bp in BPS)]
    disc = {t: [bp for bp in BPS if t in sig[bp]][0] for t in keep}
    dx = np.array([float(E[disc[t]].loc[t, "effect_size"]) for t in keep])
    dy = np.array([np.mean([float(E[bp].loc[t, "effect_size"]) for bp in BPS if bp != disc[t]])
                   for t in keep])
    axB = fig.add_subplot(gs[0, 1])
    lim = max(np.abs(dx).max(), np.abs(dy).max()) * 1.18
    axB.plot([-lim, lim], [-lim, lim], ls="--", lw=0.8, color="0.5")
    axB.axhline(0, lw=0.7, color="0.65"); axB.axvline(0, lw=0.7, color="0.65")
    axB.scatter(dx, dy, s=26, color="#1F4E79", edgecolor="white", linewidth=0.5, zorder=3)
    axB.set_xlim(-lim, lim); axB.set_ylim(-lim, lim)
    axB.set_xlabel("Effect in the cohort where\nit was significant")
    axB.set_ylabel("Mean effect in the\nother two cohorts")
    same = int((np.sign(dx) == np.sign(dy)).sum())
    axB.text(0.03, 0.97, f"n = {len(keep)} genera measurable\nin all three cohorts\n"
                         f"same direction: {same}/{len(keep)}",
             transform=axB.transAxes, ha="left", va="top", fontsize=6.5,
             bbox=dict(fc="white", ec="0.8", lw=0.5, pad=1.6))
    axB.spines[["top", "right"]].set_visible(False)
    axB.annotate("B", xy=(0, 1), xycoords="axes fraction", xytext=(-40, 6),
                 textcoords="offset points", fontsize=11, fontweight="bold", annotation_clip=False)

    axC = fig.add_subplot(gs[1, 1])
    for a, b in zip(np.abs(dx), np.abs(dy)):
        axC.plot([0, 1], [a, b], color="0.72", lw=0.8, zorder=1)
    axC.scatter(np.zeros(len(dx)), np.abs(dx), s=22, color=C_METS, edgecolor="white",
                linewidth=0.4, zorder=3)
    axC.scatter(np.ones(len(dy)), np.abs(dy), s=22, color=C_CTRL, edgecolor="white",
                linewidth=0.4, zorder=3)
    axC.plot([-0.16, 0.16], [np.abs(dx).mean()] * 2, color=C_METS, lw=2.2)
    axC.plot([0.84, 1.16], [np.abs(dy).mean()] * 2, color=C_CTRL, lw=2.2)
    axC.set_xticks([0, 1])
    axC.set_xticklabels(["Discovering\ncohort", "Other two\ncohorts"], fontsize=7)
    axC.set_xlim(-0.4, 1.4)
    axC.set_ylabel("|log fold-change|")
    drop = 100 * (1 - np.abs(dy).mean() / np.abs(dx).mean())
    axC.text(0.5, 0.97, f"mean {np.abs(dx).mean():.2f} → {np.abs(dy).mean():.2f}\n({drop:.0f}% smaller)",
             transform=axC.transAxes, ha="center", va="top", fontsize=6.6,
             bbox=dict(fc="white", ec="0.8", lw=0.5, pad=1.6))
    axC.spines[["top", "right"]].set_visible(False)
    axC.annotate("C", xy=(0, 1), xycoords="axes fraction", xytext=(-40, 6),
                 textcoords="offset points", fontsize=11, fontweight="bold", annotation_clip=False)

    os.makedirs(f"{BASE}/supplementary/optional", exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(f"{BASE}/supplementary/optional/Figure4_single_study_nonreplication.{ext}", bbox_inches="tight")
    plt.close(fig)
    print(f"Figure 4 (non-replication) written | {len(order)} genera, {len(keep)} measurable in all three, "
          f"mean |effect| {np.abs(dx).mean():.3f} -> {np.abs(dy).mean():.3f} ({drop:.0f}% drop), same direction {same}/{len(keep)}")

figure_nonreplication()
