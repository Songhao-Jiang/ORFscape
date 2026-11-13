# -*- coding: utf-8 -*-
"""TP P13 第 2/3 件：质控分布图 + 新基因事件统计。

输入：`P13_2_qc_per_peptide.tsv`（由 tp_p13_tables.py 生成）+ P7 表
输出：
  P13_2_qc_summary.tsv        质控中位数/IQR 汇总
  P13_2_qc_<metric>.png/pdf   三张分布图（cos / b-y 覆盖度 / Raw score）
  P13_3_events.tsv            事件类型 × 酶切 × ORF/肽段/谱图
  P13_3_events_per_peptide.tsv 每肽事件归属
"""
import argparse
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


from orfscape.core import novgene as G
from orfscape.core.novel import iln


def setup_font():
    """Nature 风格：拉丁 Arial + 中文 SimHei 回退。

    必须放在 import orfscape.* **之后**——`report.py` 在导入时会
    `plt.rcParams.update(...)` 把字体重置回默认，先设会被冲掉。
    """
    from matplotlib import font_manager
    fams = []
    for f in ("C:/Windows/Fonts/arial.ttf", "C:/Windows/Fonts/simhei.ttf"):
        if os.path.exists(f):
            font_manager.fontManager.addfont(f)
            fams.append(font_manager.FontProperties(fname=f).get_name())
    if fams:
        # matplotlib ≥3.6 才支持逐字回退；必须把多字体写进 font.family 本身
        plt.rcParams["font.family"] = fams + ["DejaVu Sans"]
        plt.rcParams["font.sans-serif"] = fams + ["DejaVu Sans"]
    plt.rcParams["axes.unicode_minus"] = False
    plt.rcParams["pdf.fonttype"] = 42          # 矢量字体，投稿友好
    plt.rcParams["ps.fonttype"] = 42
    # 统计量记号用 mathtext 斜体 P，字体跟随 Arial
    plt.rcParams["mathtext.fontset"] = "custom"
    plt.rcParams["mathtext.rm"] = plt.rcParams["font.sans-serif"][0]
    plt.rcParams["mathtext.it"] = plt.rcParams["font.sans-serif"][0] + ":italic"
    plt.rcParams["mathtext.bf"] = plt.rcParams["font.sans-serif"][0] + ":bold"


setup_font()

ENZ = {"TL": "Trypsin", "LL": "LysargiNase"}
ENZ_EN = {"TL": "Trypsin", "LL": "LysargiNase", "ALL": "Combined"}
SUB_EN = {"基因间区新基因": "Intergenic novel gene",
          "已注释基因不同框新基因": "Internal alternative frame",
          "反义链表达新基因": "Antisense translation",
          "已注释基因 N 端校正": "N-terminal correction",
          "未定位": "Unassigned"}
TYPE_ZH = {
    "intergenic new coding gene": "基因间区新基因",
    "annotated gene correction": "已注释基因 N 端校正",
}
# classify_events 的内部细分键 → 中文
SUB_ZH = {
    "intergenic_only": "基因间区新基因",
    "same_strand_other_frame": "已注释基因不同框新基因",
    "antisense_overlap": "反义链表达新基因",
    "nterm_extension": "已注释基因 N 端校正",
    "unmapped": "未定位",
}


def zh_sub(s):
    return SUB_ZH.get(str(s).strip().lower(), str(s))


def stars_for(p):
    if p < 1e-4:
        return "****"
    if p < 1e-3:
        return "***"
    if p < 1e-2:
        return "**"
    if p < 5e-2:
        return "*"
    return "n.s."


def stars_test(a, b):
    """两组比较：双侧 Mann–Whitney U（Wilcoxon 秩和，含并列校正）。

    Nature 图注口径：*P<0.05, **P<0.01, ***P<0.001, ****P<0.0001，n.s. 不显著。
    分布明显偏态（cos 双峰、Raw score 长尾），不用 t 检验。
    """
    from scipy import stats as _st
    try:
        u, p = _st.mannwhitneyu(a, b, alternative="two-sided")
    except ValueError:
        return float("nan"), float("nan")
    return float(u), float(p)


def add_sig(ax, x1, x2, y, p, h=None):
    """Nature 风格显著性括号：横线 + 两竖 + 星号 + p 值。"""
    span = ax.get_ylim()[1] - ax.get_ylim()[0]
    h = h if h is not None else span * 0.035
    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=0.8, color="black",
            clip_on=False)
    ax.text((x1 + x2) / 2, y + h * 1.15, stars_for(p), ha="center", va="bottom",
            fontsize=9, fontweight="bold", color="black")
    ax.text((x1 + x2) / 2, y + h * 1.15 - span * 0.055,
            ("$P$ = %.2g" % p) if p >= 1e-4 else "$P$ < 0.0001",
            ha="center", va="top", fontsize=6, color="#333333")


def boxplot(ax, data, labels, colors, ylabel, title, p=None, ylim=None):
    """箱线 + 抖动散点 + Mann–Whitney U 显著性标注。"""
    bp = ax.boxplot(data, labels=labels, widths=0.5, patch_artist=True,
                    showfliers=False, whis=1.5,
                    boxprops=dict(lw=0.7, edgecolor="black"),
                    whiskerprops=dict(lw=0.7, color="black"),
                    capprops=dict(lw=0.7, color="black"),
                    medianprops=dict(color="black", lw=1.1))
    for i, patch in enumerate(bp["boxes"]):
        patch.set_facecolor(colors[i % len(colors)])
        patch.set_alpha(0.55)
    # 抖动散点（Nature 常见做法：把每个肽画出来）
    rng = np.random.default_rng(7)
    for i, vals in enumerate(data):
        v = np.asarray(vals, dtype=float)
        v = v[np.isfinite(v)]
        if not len(v):
            continue
        x = np.full(len(v), i + 1) + rng.uniform(-0.14, 0.14, len(v))
        ax.scatter(x, v, s=3.2, color="#4B5563", alpha=0.35, linewidths=0,
                   zorder=1)
    ax.set_ylabel(ylabel, fontsize=7.5)
    ax.set_title(title, fontsize=8, loc="left", fontweight="bold")
    ax.tick_params(labelsize=7, width=0.7, length=2.5)
    ax.tick_params(axis="x", length=0, pad=2)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_linewidth(0.7)
    if p is not None and np.isfinite(p):
        vals = [np.asarray(d, dtype=float) for d in data]
        vals = [v[np.isfinite(v)] for v in vals if np.isfinite(v).any()]
        ymax = max([v.max() for v in vals] or [1])
        ymin = min([v.min() for v in vals] or [0])
        if ylim is not None:
            top = ylim[1]
            ax.set_ylim(*ylim)
            add_sig(ax, 1, 2, top - (top - ylim[0]) * 0.14, p)
        else:
            pad = (ymax - ymin) * 0.10
            ax.set_ylim(ymin - (ymax - ymin) * 0.05, ymax + pad * 2.0)
            add_sig(ax, 1, 2, ymax + pad * 0.55, p)
    elif ylim is not None:
        ax.set_ylim(*ylim)


def main(argv=None):
    ap = argparse.ArgumentParser(description="TP P13 质控图与事件统计")
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--res", required=True, help="P13_results 目录")
    args = ap.parse_args(argv)

    root = os.path.realpath(args.workdir)
    res = os.path.realpath(args.res)
    qc = pd.read_csv(os.path.join(res, "P13_2_qc_per_peptide.tsv"),
                     sep="\t", encoding="utf-8-sig")
    qc["enzyme"] = qc["enzyme"].fillna("")
    # 用 isoPSV + 补算合成的余弦，并只保留有值的肽 → 四格指标同一套肽（数字对齐）
    if "cosine_final" in qc.columns:
        qc = qc[qc["cosine_final"].notna()].copy()
        qc["cosine"] = qc["cosine_final"]

    # ---------- 质控汇总
    rows = []
    for cat in ("Filtered", "Obtained"):
        for enz in ("Trypsin", "LysargiNase", "Combined"):
            s = qc[qc["category"] == cat]
            if enz != "Combined":
                s = s[s["enzyme"] == enz]
            if not len(s):
                continue
            def q(col):
                x = pd.to_numeric(s[col], errors="coerce").dropna()
                if not len(x):
                    return (float("nan"),) * 3
                return (float(np.median(x)),
                        float(np.percentile(x, 25)),
                        float(np.percentile(x, 75)))
            c = q("cosine")
            v = q("ion_coverage")
            r = q("raw_score")
            sp = q("n_spectra")
            rows.append(dict(category=cat, group=enz, n=len(s),
                             cosine_median=c[0], cosine_Q1=c[1], cosine_Q3=c[2],
                             ion_coverage_median=v[0], ion_coverage_Q1=v[1],
                             ion_coverage_Q3=v[2],
                             raw_score_median=r[0], raw_score_Q1=r[1],
                             raw_score_Q3=r[2],
                             spectra_median=sp[0], spectra_Q1=sp[1], spectra_Q3=sp[2]))
    summ = pd.DataFrame(rows)
    summ.to_csv(os.path.join(res, "P13_2_qc_summary.tsv"), sep="\t",
                index=False, encoding="utf-8-sig")

    # ---------- 质控图：Filtered vs Obtained（4 指标 × 3 组）
    specs = [("cosine", "Cosine similarity (AlphaPeptDeep)", "Cosine similarity", (0, 1.0)),
             ("ion_coverage", "b/y ion coverage (±20 ppm)", "Ion coverage", (0, 1.0)),
             ("raw_score", "pFind Raw score", "Raw score", None),
             ("n_spectra", "Spectra per peptide", "PSMs per peptide", None)]
    groups = [("Trypsin", "A  Trypsin"), ("LysargiNase", "B  LysargiNase"),
              ("ALL", "C  Combined")]
    cats = ("Filtered", "Obtained")
    cols = {"Filtered": "#C0504D", "Obtained": "#2F6FBA"}
    stat_rows = []
    fig, axes = plt.subplots(len(specs), len(groups),
                             figsize=(7.2, 1.85 * len(specs)))
    for r, (key, title, ylab, ylim) in enumerate(specs):
        for cix, (enz, panel_title) in enumerate(groups):
            ax = axes[r][cix]
            data, labels, colors = [], [], []
            for cat in cats:
                s = qc[qc["category"] == cat]
                if enz != "ALL":
                    s = s[s["enzyme"] == enz]
                vals = pd.to_numeric(s[key], errors="coerce").dropna().values
                data.append(vals if len(vals) else np.array([np.nan]))
                labels.append(("%s" + chr(10) + "n=%d") % (cat, len(vals)))
                colors.append(cols[cat])
            a, b = data
            if len(a) > 2 and len(b) > 2:
                u, p = stars_test(a, b)
            else:
                u, p = float("nan"), float("nan")
            stat_rows.append(dict(metric=key, group=enz,
                                  n_filtered=len(a), n_obtained=len(b),
                                  median_filtered=float(np.nanmedian(a)) if len(a) else float("nan"),
                                  median_obtained=float(np.nanmedian(b)) if len(b) else float("nan"),
                                  U=u, p=p, stars=stars_for(p)))
            boxplot(ax, data, labels, colors, ylab if cix == 0 else "",
                    "", p=p, ylim=ylim)
            if key == "n_spectra":      # 谱数跨三个数量级，用对数轴
                ax.set_yscale("log")
                hi = float(np.nanmax(np.concatenate([a, b]))) if len(a) + len(b) else 10
                ax.set_ylim(0.8, max(10.0, hi * 1.8))
            if r == 0:
                ax.set_title(panel_title, fontsize=8, loc="left", fontweight="bold")
    fig.suptitle("Novel peptide quality control: Filtered vs Obtained", fontsize=9, y=0.995)
    fig.text(0.5, 0.004,
             "Two-sided Mann–Whitney U test; *P<0.05, **P<0.01, ***P<0.001, "
             "****P<0.0001; n.s., not significant. Box, median and IQR; "
             "whiskers, 1.5×IQR; dots, single peptides. Cosine is "
             "AlphaPeptDeep (isoPSV values, supplemented by in-house calculation "
             "for peptides isoPSV could not score); all panels use the same "
             "peptide set; 21 Filtered peptides without any prediction are excluded.",
             ha="center", va="bottom", fontsize=5.6, color="#333333")
    fig.tight_layout(rect=(0, 0.03, 1, 0.97))
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(res, "P13_2_qc_panels.%s" % ext),
                    dpi=600 if ext == "png" else 300, bbox_inches="tight")
    plt.close(fig)
    pd.DataFrame(stat_rows).to_csv(os.path.join(res, "P13_2_qc_stats.tsv"),
                                   sep="	", index=False, encoding="utf-8-sig")

    # ---------- 事件分类（按 KEEP 行）
    v = pd.read_csv(os.path.join(root, "06b_validate_research", "novel_validated.tsv"),
                    sep="\t", dtype=str, encoding="utf-8-sig")
    v["pep"] = v["peptide"].map(iln)
    keep = v[v["final_keep"].astype(str).str.startswith("KEEP")].copy()
    ev = G.classify_events(root, keep.to_dict("records"), log=print)
    ev = pd.DataFrame(ev)
    ev["pep"] = ev["peptide"].map(iln)
    ev["sub"] = ev["event_sub"].map(zh_sub)
    ev["event_label"] = ev["sub"].map(SUB_EN)
    # 输出全英文：去掉中文 sub 列（内部逻辑已用完）
    ev.drop(columns=["sub"]).to_csv(
        os.path.join(res, "P13_3_events_per_peptide.tsv"), sep="\t",
        index=False, encoding="utf-8-sig")

    # 每肽谱数（分酶）
    psm = defaultdict(lambda: defaultdict(int))
    for sub in ("04_search_disc", "05b_research"):
        base = os.path.join(root, sub, "TP")
        if not os.path.isdir(base):
            continue
        from orfscape.core.novel import load_spectra
        for t in sorted(os.listdir(base)):
            f = os.path.join(base, t, "pFind-Filtered.spectra")
            if not os.path.exists(f):
                continue
            enz = "LL" if "_LL_" in t else "TL"
            for r in load_spectra(f):
                p = iln(r.get("Sequence") or "")
                if p:
                    psm[p][enz] += 1

    rows = []
    for enz in ("TL", "LL", "ALL"):
        for sub_zh in ("基因间区新基因", "已注释基因不同框新基因",
                       "反义链表达新基因", "已注释基因 N 端校正", "未定位"):
            s = ev[ev["sub"] == sub_zh]
            if enz != "ALL":
                s = s[s["pep"].map(lambda p: psm[p].get(enz, 0) > 0)]
            if not len(s):
                continue
            n = sum(sum(psm[p].values()) if enz == "ALL" else psm[p].get(enz, 0)
                    for p in s["pep"])
            rows.append(dict(Enzyme=ENZ_EN.get(enz, "Combined"),
                             Event_type=SUB_EN[sub_zh],
                             ORFs=s["event_orf"].nunique(),
                             Peptides=s["pep"].nunique(), PSMs=int(n)))
    t3 = pd.DataFrame(rows)
    t3.to_csv(os.path.join(res, "P13_3_events.tsv"), sep="\t",
              index=False, encoding="utf-8-sig")
    print("\n== 事件统计 ==")
    print(t3.to_string(index=False))
    print("\n== 质控汇总 ==")
    print(summ.to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
