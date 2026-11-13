# -*- coding: utf-8 -*-
"""R9 Nature 产出器：主图 + 附表 + Methods + Results（自动填充全部数字）。

图规格：matplotlib，Arial 7pt，89mm 单栏 / 183mm 双栏，PDF+SVG+PNG 300dpi。
表规格：xlsx 多 sheet，无冗饰。
"""
import csv
import datetime
import json
import os
from collections import Counter, defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .. import __version__
from .inputs import project_species, safe_segment

# 物种配色：按项目物种顺序循环，不写死项目专属代码
_PALETTE = ("#2F6FB0", "#E07B39", "#4C9A2A", "#8E5AA8", "#C0392B", "#2A9D8F")


def _sp_color(i):
    return _PALETTE[i % len(_PALETTE)]


def _novels(root):
    """{物种: 新肽表}；物种取自项目本身，不写死项目专属代码。"""
    out = {}
    for sp in project_species(root):
        seg = safe_segment(sp)
        out[seg] = load_tsv(os.path.join(root, "05_novel", seg, "novel_peptides.tsv"))
    return out

# Nature 风格
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans"],
    "font.size": 7,
    "axes.linewidth": 0.5,
    "axes.labelsize": 7,
    "xtick.labelsize": 6,
    "ytick.labelsize": 6,
    "legend.fontsize": 6,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
})
MM = 1 / 25.4  # mm → inch
COL1 = 89 * MM   # 单栏
COL2 = 183 * MM  # 双栏

# 色盲安全
COLORS = ["#0173B2", "#DE8F05", "#029E73", "#D55E00", "#CC78BC",
          "#CA9161", "#FBAFE4", "#949494", "#ECE133", "#56B4E9"]

# ---------------- i18n（默认英文；--lang zh 切中文）----------------
from .i18n import get_lang

ZH_LABELS = {
    "Annotated DB search": "注释库搜库", "PSMs": "PSM 数",
    "PSM distribution": "PSM 分布", "Database comparison": "数据库对比",
    "Avg PSMs / task": "平均 PSM/任务", "Injection time": "注入时间",
    "MS2 IT saturation (%)": "MS2 注入时间触顶率 (%)", "Raw files": "原始文件数",
    "Cleavage specificity": "酶切特异性", "Fully cleaved (%)": "全切率 (%)",
    "Trypsin+LysC": "胰酶+LysC", "LysargiNase+LysN": "LysargiNase+LysN",
    "AGC fill": "AGC 填充", "Median MS2 AGC fill": "MS2 AGC 填充中位数",
    "E2 PSMs": "E2 PSM 数", "log10(count)": "log10(计数)",
    "By species": "按物种", "Novel peptides": "新肽段", "Source layer": "来源层",
    "Tier1 mirror peptides": "Tier1 镜像肽段", "Peptide length (aa)": "肽段长度 (aa)",
    "Count": "计数", "PSM support by tier": "各 Tier 的 PSM 支持",
    "PSM count (capped at 50)": "PSM 数（上限 50）",
    "intergenic_6frame": "基因间六框", "alt_start_extension": "替代起始延伸",
    "unknown": "未知",
    "a  Annotated DB search": "a  注释库搜库",
    "b  PSM distribution": "b  PSM 分布",
    "c  Database comparison": "c  数据库对比",
    "a  Injection time": "a  注入时间",
    "b  Cleavage specificity": "b  酶切特异性",
    "c  AGC fill": "c  AGC 填充",
    "a  Discovery funnel": "a  发现漏斗",
    "b  By species": "b  按物种", "c  Source layer": "c  来源层",
    "a  Tier1 mirror peptides": "a  Tier1 镜像肽段",
    "b  PSM support by tier": "b  各 Tier 的 PSM 支持",
    # 表头
    "Species": "物种", "Enzyme": "酶", "Tasks": "任务数", "Total_PSM": "总 PSM",
    "Task": "任务", "PSM": "PSM", "Peptides": "肽段",
    "Full_cleavage_%": "全切率_%", "MC0_%": "零漏切_%",
    "Parameter": "参数", "Value": "值", "Tier": "Tier", "Layer": "层",
    "Locus": "基因座", "Enzymes": "酶",
}

def _t(en):
    """图表/表头标签翻译：lang=zh 返回中文，否则英文原文。"""
    return ZH_LABELS.get(en, en) if get_lang() == "zh" else en

def _apply_zh_font():
    """zh 模式切换中文字体（微软雅黑），避免 CJK 缺字方框。"""
    if get_lang() == "zh":
        plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]



def load_tsv(path):
    if not os.path.exists(path):
        return []
    with open(path, newline="", encoding="utf-8-sig") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def build_all(root, out_dir):
    """生成全部 Nature 产出。返回产出文件列表。"""
    os.makedirs(out_dir, exist_ok=True)
    _apply_zh_font()
    files = []
    files += fig1_workflow(root, out_dir)
    files += fig2_qc(root, out_dir)
    files += fig3_novel_funnel(root, out_dir)
    files += fig4_mirror(root, out_dir)
    files += tables(root, out_dir)
    files += methods_md(root, out_dir)
    return files


def _save(fig, out_dir, name):
    paths = []
    for ext in ("pdf", "svg", "png"):
        p = os.path.join(out_dir, f"{name}.{ext}")
        fig.savefig(p, format=ext)
        paths.append(p)
    plt.close(fig)
    return paths


def _save_tsv(out_dir, name, header, rows):
    """与 _save 同款的落盘写法（表：制表符分隔、utf-8-sig）。"""
    p = os.path.join(out_dir, f"{name}.tsv")
    with open(p, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(header)
        w.writerows(rows)
    return p


# ─── 补充图 EVID：个体级证据（理论 × 实测）──────────────────────
MIN_SCORE = 0.28
MIN_COS = 0.82
_VALID_DIRS = ("06_validate", "06b_validate_research")


def _load_evidence_points(root):
    """读 <validate_dir>/novel_validated.tsv（scandir 白名单）→ 点位列表。"""
    name = os.environ.get("ORFSCAPE_VALIDATE_DIR", "06_validate")
    if name not in _VALID_DIRS:
        raise ValueError(f"unsafe validate dir: {name}")
    d = os.path.join(str(root), name)
    pts = []
    if not os.path.isdir(d):
        return pts
    with os.scandir(d) as it:
        for e in it:
            if e.name != "novel_validated.tsv" or not e.is_file():
                continue
            with open(e.path, encoding="utf-8-sig", newline="") as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    try:
                        cos = float(r.get("v_best_cos") or "")
                        score = float(r.get("v_best_score") or "")
                    except (TypeError, ValueError):
                        continue
                    keep = str(r.get("final_keep", "")).startswith("KEEP")
                    pts.append(dict(peptide=r["peptide"],
                                    protein=r.get("protein", ""),
                                    species=r.get("species", ""),
                                    tier=r.get("tier_final") or r.get("tier", ""),
                                    cos=cos, score=score, keep=keep,
                                    reason=r.get("keep_reason", "") if keep else "",
                                    psm=r.get("psm", ""), is_new=r.get("is_new", ""),
                                    doubt=r.get("v_doubt", "")))
    return pts


def fig_evid1_theory_vs_observed(root, out_dir):
    """每点一条新肽：横轴预测-实测余弦（理论），纵轴 isoPSV cand1 打分（实测）。"""
    pts = _load_evidence_points(root)
    if not pts:
        return []
    out_hi = [p for p in pts if p["keep"] and p["cos"] < MIN_COS
              and p["score"] >= MIN_SCORE]
    out_lo = [p for p in pts if p["keep"] and p["cos"] >= MIN_COS
              and p["score"] < MIN_SCORE]
    colors = {"isoPSV+mirror": ("#2F6FB0", "isoPSV + mirror"),
              "mirror_only": ("#E07B39", "mirror only (isoPSV failed)"),
              "isoPSV": ("#3E9B62", "isoPSV (Tier3, 3-gate)"),
              "": ("#B8B8B8", "dropped")}
    fig, ax = plt.subplots(figsize=(3.4, 3.0))
    for key, (color, label) in colors.items():
        sub = [p for p in pts if p["reason"] == key]
        if not sub:
            continue
        ax.scatter([p["cos"] for p in sub], [p["score"] for p in sub],
                   s=5, c=color, alpha=0.55, linewidths=0,
                   label=f"{label} (n={len(sub)})", rasterized=True)
    ax.axvline(MIN_COS, color="#666666", lw=0.6, ls="--")
    ax.axhline(MIN_SCORE, color="#666666", lw=0.6, ls="--")
    ax.set_xlim(-0.02, 1.02)
    ax.set_xlabel("Predicted-observed spectral cosine (theoretical)")
    ax.set_ylabel("isoPSV cand1 score (observed)")
    ax.set_title("Per-peptide evidence (one dot = one novel peptide)",
                 fontsize=8, pad=4)
    ax.text(0.98, 0.97, "high theory / high obs", transform=ax.transAxes,
            fontsize=5.2, ha="right", va="top", color="#555555")
    ax.text(0.02, 0.97, "low theory / high obs", transform=ax.transAxes,
            fontsize=5.2, ha="left", va="top", color="#B03A2E")
    ax.text(0.02, 0.91, f"n={len(out_hi)}", transform=ax.transAxes,
            fontsize=6.5, ha="left", va="top", color="#B03A2E", weight="bold")
    ax.text(0.98, 0.03, "high theory / low obs", transform=ax.transAxes,
            fontsize=5.2, ha="right", va="bottom", color="#555555")
    ax.text(0.02, 0.03, "low / low", transform=ax.transAxes,
            fontsize=5.2, ha="left", va="bottom", color="#555555")
    ax.legend(loc="lower right", fontsize=5, frameon=False,
              handletextpad=0.4, borderpad=0.2, bbox_to_anchor=(1.0, 0.10))
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    rows = []
    for p in out_hi:
        rows.append(["low_theory_high_obs", p["peptide"], p["protein"],
                     p["species"], p["tier"], p["cos"], p["score"],
                     p["reason"], p["psm"], p["is_new"], p["doubt"]])
    for p in out_lo:
        rows.append(["high_theory_low_obs", p["peptide"], p["protein"],
                     p["species"], p["tier"], p["cos"], p["score"],
                     p["reason"], p["psm"], p["is_new"], p["doubt"]])
    files = _save(fig, out_dir, "Fig_EVID1_theory_vs_observed")
    files.append(_save_tsv(
        out_dir, "EVID1_outliers",
        ["class", "peptide", "protein", "species", "tier", "cos", "score",
         "keep_reason", "psm", "is_new", "doubt"], rows))
    return files


def fig_evid2_kept_vs_dropped(root, out_dir):
    pts = _load_evidence_points(root)
    if not pts:
        return []
    fig, axes = plt.subplots(1, 2, figsize=(4.6, 1.9))
    for ax, key, xlab in ((axes[0], "cos", "Predicted-observed cosine"),
                          (axes[1], "score", "isoPSV cand1 score")):
        keep = [p[key] for p in pts if p["keep"]]
        drop = [p[key] for p in pts if not p["keep"]]
        bins = [i / 40.0 for i in range(41)] if key == "cos" else 20
        ax.hist([keep, drop], bins=bins,
                label=[f"kept (n={len(keep)})", f"dropped (n={len(drop)})"],
                color=["#2F6FB0", "#B8B8B8"], histtype="step", lw=0.9)
        ax.set_xlabel(xlab)
        ax.set_ylabel("Novel peptides")
        ax.legend(fontsize=5, frameon=False)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    fig.tight_layout()
    return _save(fig, out_dir, "Fig_EVID2_kept_vs_dropped")


def _count_faa(path):
    try:
        with open(path, encoding="utf-8", errors="replace") as fh:
            return sum(1 for ln in fh if ln.startswith(">"))
    except OSError:
        return 0


def _novel_yield(root, tag):
    """对 E2 的某一库（aug/6f）跑四滤，返回 {(sp): 唯一新肽数}。"""
    from .novel import (digest_krc, digest_krn, isobaric_variants, iln,
                        load_spectra, read_faa)
    out = {}
    for sp in project_species(root):
        db = os.path.join(root, "02_db", sp)
        recs = read_faa(os.path.join(db, "db_annot.faa"))
        l0 = [(h, s) for h, s in recs if h.startswith("L0_")]
        known = set()
        for _, s in l0:
            known.update(digest_krc(s))
            known.update(digest_krn(s))
        big = "|".join(s for _, s in l0)
        e1 = set()
        ap = os.path.join(root, "04_search_annot", sp)
        if os.path.isdir(ap):
            for t in os.listdir(ap):
                p = os.path.join(ap, t, "pFind-Filtered.spectra")
                if os.path.exists(p):
                    for r in load_spectra(p):
                        e1.add(iln(r.get("Sequence", "")))
        pep = set()
        base = os.path.join(root, "04_search_disc", sp)
        if os.path.isdir(base):
            for t in os.listdir(base):
                if not t.endswith("__" + tag):
                    continue
                p = os.path.join(base, t, "pFind-Filtered.spectra")
                if not os.path.exists(p):
                    continue
                for r in load_spectra(p):
                    pr = r.get("Proteins") or ""
                    q = iln(r.get("Sequence", ""))
                    if not q or not ("L2_" in pr or "L3_" in pr):
                        continue
                    if q in known or q in e1 or q in big:
                        continue
                    if any(v in known or v in e1
                           for v in isobaric_variants(q) if v != q):
                        continue
                    pep.add(q)
        out[sp] = len(pep)
    return out


def fig_evid3_layer_gain(root, out_dir):
    """库层边际增益：三档库（annot / +Prodigal+N端 / +六框）的条目数与新肽。"""
    sizes, novel = {}, {}
    for sp in project_species(root):
        db = os.path.join(root, "02_db", sp)
        n_annot = _count_faa(os.path.join(db, "db_annot.faa"))
        n_l1 = _count_faa(os.path.join(db, "L1_prodigal_filt.faa"))
        n_l3 = _count_faa(os.path.join(db, "L3_altstart.faa"))
        sizes[sp] = (n_annot, n_annot + n_l1 + n_l3,
                     _count_faa(os.path.join(db, "db_6f.faa")))
    y_aug = _novel_yield(root, "aug")
    y_6f = _novel_yield(root, "6f")
    layers = ("annotated", "+Prodigal/N-term", "+six-frame")
    fig, ax = plt.subplots(figsize=(3.6, 2.5))
    width = 0.36
    for k, sp in enumerate(project_species(root)):
        color = _sp_color(k)
        ys = (0, y_aug.get(sp, 0), y_6f.get(sp, 0))
        xs = [i + (k - 0.5) * width for i in range(3)]
        ax.bar(xs, ys, width=width, color=color, label=sp, edgecolor="white",
               linewidth=0.5)
        for i, (x, y) in enumerate(zip(xs, ys)):
            ax.text(x, y + max(40, max(ys) * 0.03), f"{y:,}", ha="center",
                    va="bottom", fontsize=5, color="#444444")
    for i, nm in enumerate(layers):
        ax.text(i, -max(max(y_6f.values()), 1) * 0.16,
                "n=%s" % "/".join(f"{sizes[sp][i]:,}" for sp in ("Mab", "Mch")),
                ha="center", va="top", fontsize=4.6, color="#777777")
    ax.set_xticks(range(3))
    ax.set_xticklabels(layers, fontsize=6)
    ax.set_ylabel("Unique novel peptides (after 4 filters)")
    ax.set_title("Marginal gain of each database layer", fontsize=8, pad=16)
    ax.legend(fontsize=5.5, frameon=False, loc="upper left")
    ax.set_ylim(0, max(max(y_6f.values()), 1) * 1.18)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    return _save(fig, out_dir, "Fig_EVID3_layer_gain")


def fig_evid4_degraded_branch(root, out_dir):
    """无镜像数据时的降级分支：单酶 → 仅 Tier3 → isoPSV 三门。"""
    from matplotlib.patches import FancyArrowPatch, Rectangle

    fig, ax = plt.subplots(figsize=(5.2, 3.4))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis("off")

    def box(x, y, w, h, text, fc="#EEF3F8", ec="#2F6FB0", fs=6.2):
        ax.add_patch(Rectangle((x, y), w, h, facecolor=fc, edgecolor=ec,
                               linewidth=0.7, zorder=2))
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center",
                fontsize=fs, zorder=3, linespacing=1.35)

    def arrow(x1, y1, x2, y2, style="-|>", color="#555555", ls="-"):
        ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle=style,
                                     mutation_scale=7, linewidth=0.6,
                                     color=color, linestyle=ls, zorder=1))

    box(3.4, 8.9, 3.2, 0.8, "raw + species\n(auto QC, FAIMS split)", "#F2F2F2", "#888888")
    box(3.4, 7.6, 3.2, 0.8, "4 filters -> GSSP\n(locus-resolved novel peptides)")
    arrow(5.0, 8.9, 5.0, 8.4)
    box(3.4, 6.3, 3.2, 0.8, "P6b candidate-ORF re-search\n(6k DB, same MGFs)")
    arrow(5.0, 7.6, 5.0, 7.1)

    # 分支
    arrow(4.3, 6.3, 2.5, 5.5)
    arrow(5.7, 6.3, 7.5, 5.5)
    box(0.6, 4.5, 3.8, 1.0,
        "SINGLE ENZYME (no mirror)\nall peptides -> Tier 3", "#FDF0E6", "#E07B39")
    box(5.6, 4.5, 3.8, 1.0,
        "MIRROR PAIR (TL + LL)\nTier 1 / Tier 2 / Tier 3", "#EAF4EE", "#3E9B62")

    box(0.6, 3.0, 3.8, 0.9,
        "isoPSV 3 gates\nisoPSV x quality x similarity", "#FDF0E6", "#E07B39")
    box(5.6, 3.0, 3.8, 0.9,
        "Tier 1/2: keep (mirror evidence)\nTier 3: isoPSV 3 gates", "#EAF4EE", "#3E9B62")
    arrow(2.5, 4.5, 2.5, 3.9)
    arrow(7.5, 4.5, 7.5, 3.9)

    box(3.4, 1.8, 3.2, 0.8, "KEEP\n(reason recorded per peptide)", "#EEF3F8", "#2F6FB0")
    arrow(2.5, 3.0, 4.4, 2.6)
    arrow(7.5, 3.0, 5.6, 2.6)
    box(3.4, 0.5, 3.2, 0.8, "Event classification\n-> structure / function / homology")
    arrow(5.0, 1.8, 5.0, 1.3)

    ax.text(2.5, 5.82, "degraded path", fontsize=5.4, color="#E07B39",
            ha="center", style="italic")
    ax.text(7.5, 5.82, "standard path", fontsize=5.4, color="#3E9B62",
            ha="center", style="italic")
    fig.suptitle("Degraded branch when no mirror-enzyme data exist",
                 fontsize=8, y=0.995)
    fig.subplots_adjust(left=0.01, right=0.99, top=0.93, bottom=0.01)
    return _save(fig, out_dir, "Fig_EVID4_degraded_branch")


# ─── Fig 1: 工作流+鉴定全景 ───
def fig1_workflow(root, out_dir):
    e1 = load_tsv(os.path.join(root, "04_search_annot", "e1_report.tsv"))
    e2 = load_tsv(os.path.join(root, "04_search_disc", "e2_report.tsv"))
    if not e1:
        return []
    fig, axes = plt.subplots(1, 3, figsize=(COL2, COL2 * 0.22))

    # a: 分组 PSM
    ax = axes[0]
    groups = defaultdict(int)
    for r in e1:
        groups[(r["species"], r["enzyme"])] += int(r["n_psm"])
    labels = [f"{k[0]}\n{k[1]}" for k in sorted(groups)]
    vals = [groups[k] for k in sorted(groups)]
    ax.bar(range(len(vals)), vals, color=COLORS[:len(vals)], edgecolor="black", linewidth=0.3)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, fontsize=5)
    ax.set_ylabel(_t("PSMs"))
    ax.set_title("a  Annotated DB search", fontsize=7, fontweight="bold", loc="left")

    # b: 每任务 PSM
    ax = axes[1]
    psm_vals = [int(r["n_psm"]) for r in e1 if r["status"] == "OK"]
    ax.hist(psm_vals, bins=30, color=COLORS[0], edgecolor="black", linewidth=0.3)
    ax.set_xlabel(_t("PSMs per task"))
    ax.set_ylabel(_t("Count"))
    ax.set_title("b  PSM distribution", fontsize=7, fontweight="bold", loc="left")

    # c: 库大小 vs PSM
    ax = axes[2]
    if e2:
        dbs = defaultdict(lambda: [0, 0])
        for r in e2:
            dbs[r["db"]][0] += 1
            dbs[r["db"]][1] += int(r["n_psm"])
        db_names = sorted(dbs.keys())
        aug_psm = [dbs[d][1] / dbs[d][0] for d in db_names]
        ax.bar(range(len(db_names)), aug_psm, color=COLORS[2], edgecolor="black", linewidth=0.3)
        ax.set_xticks(range(len(db_names)))
        ax.set_xticklabels([f"db_{d}" for d in db_names], fontsize=6)
        ax.set_ylabel(_t("Avg PSMs / task"))
        ax.set_title("c  Database comparison", fontsize=7, fontweight="bold", loc="left")
    return _save(fig, out_dir, "Fig1_workflow_identification")


# ─── Fig 2: QC ───
def fig2_qc(root, out_dir):
    qc = load_tsv(os.path.join(root, "04a_qc", "qc_per_task.tsv"))
    raw_info = load_tsv(os.path.join(root, "00_input", "raw_info.tsv"))
    if not qc:
        return []
    fig, axes = plt.subplots(1, 3, figsize=(COL2, COL2 * 0.22))

    # a: IT 触顶率
    ax = axes[0]
    sat_vals = [float(r["it_sat_pct"]) for r in raw_info if r.get("it_sat_pct")]
    ax.hist(sat_vals, bins=20, color=COLORS[1], edgecolor="black", linewidth=0.3)
    ax.axvline(95, color="red", linestyle="--", linewidth=0.5)
    ax.set_xlabel(_t("MS2 IT saturation (%)"))
    ax.set_ylabel(_t("Raw files"))
    ax.set_title("a  Injection time", fontsize=7, fontweight="bold", loc="left")

    # b: 酶切特异性
    ax = axes[1]
    tl = [float(r["pct_full_cleavage"]) for r in qc if "_TL_" in r["task"]]
    ll = [float(r["pct_full_cleavage"]) for r in qc if "_LL_" in r["task"]]
    bp = ax.boxplot([tl, ll], labels=[_t("Trypsin+LysC"), "LysargiNase+LysN"], patch_artist=True)
    for patch, color in zip(bp["boxes"], [COLORS[0], COLORS[3]]):
        patch.set_facecolor(color)
        patch.set_alpha(0.6)
    ax.set_ylabel(_t("Fully cleaved (%)"))
    ax.set_title("b  Cleavage specificity", fontsize=7, fontweight="bold", loc="left")

    # c: AGC Fill
    ax = axes[2]
    fill_vals = [float(r["med_fill_ms2"]) for r in raw_info if r.get("med_fill_ms2")]
    ax.hist(fill_vals, bins=20, color=COLORS[4], edgecolor="black", linewidth=0.3)
    ax.set_xlabel(_t("Median MS2 AGC fill"))
    ax.set_ylabel(_t("Raw files"))
    ax.set_title("c  AGC fill", fontsize=7, fontweight="bold", loc="left")
    return _save(fig, out_dir, "Fig2_QC")


# ─── Fig 3: 新肽漏斗 ───
def fig3_novel_funnel(root, out_dir):
    novels = _novels(root)
    if not any(novels.values()):
        return []
    fig, axes = plt.subplots(1, 3, figsize=(COL2, COL2 * 0.22))

    # a: 漏斗
    ax = axes[0]
    e2 = load_tsv(os.path.join(root, "04_search_disc", "e2_report.tsv"))
    total_psm = sum(int(r["n_psm"]) for r in e2) if e2 else 0
    total_novel = sum(len(v) for v in novels.values())
    stages = ["E2 PSMs", "Novel\npeptides", "Tier1\nmirror", "Tier2\ndual-enzyme"]
    t1 = sum(1 for v in novels.values() for r in v if r.get("tier") == "Tier1")
    t2 = sum(1 for v in novels.values() for r in v if r.get("tier") == "Tier2")
    vals = [total_psm, total_novel, t1, t2]
    y_pos = range(len(stages))
    ax.barh(y_pos, np.log10([max(v, 1) for v in vals]), color=COLORS[:4],
            edgecolor="black", linewidth=0.3)
    ax.set_yticks(y_pos)
    _nl = chr(10)
    ax.set_yticklabels([_nl.join(_t(part) for part in s.split(_nl)) for s in stages],
                       fontsize=6)
    ax.set_xlabel(_t("log10(count)"))
    for i, v in enumerate(vals):
        ax.text(np.log10(max(v, 1)) + 0.1, i, f"{v:,}", va="center", fontsize=5)
    ax.set_title("a  Discovery funnel", fontsize=7, fontweight="bold", loc="left")

    # b: 按物种
    ax = axes[1]
    species = list(novels.keys())
    tier_counts = {}
    for sp in species:
        tc = Counter(r.get("tier", "Tier3") for r in novels[sp])
        tier_counts[sp] = tc
    tier_names = ["Tier1", "Tier2", "Tier3"]
    x = np.arange(len(species))
    w = 0.25
    for i, tn in enumerate(tier_names):
        vals = [tier_counts[sp].get(tn, 0) for sp in species]
        ax.bar(x + i * w, vals, w, label=tn, color=COLORS[i], edgecolor="black", linewidth=0.3)
    ax.set_xticks(x + w)
    ax.set_xticklabels(species)
    ax.set_ylabel(_t("Novel peptides"))
    ax.legend()
    ax.set_title("b  By species", fontsize=7, fontweight="bold", loc="left")

    # c: 来源层
    ax = axes[2]
    src = Counter()
    for sp in species:
        for r in novels[sp]:
            src[r.get("source", "unknown")] += 1
    ax.pie(src.values(), labels=[s.replace("_", "\n") for s in src.keys()],
           colors=COLORS[:len(src)], autopct="%1.0f%%", textprops={"fontsize": 5})
    ax.set_title("c  Source layer", fontsize=7, fontweight="bold", loc="left")
    return _save(fig, out_dir, "Fig3_novel_funnel")


# ─── Fig 4: 镜像验证 ───
def fig4_mirror(root, out_dir):
    novels = _novels(root)
    if not any(novels.values()):
        return []
    fig, axes = plt.subplots(1, 2, figsize=(COL2 * 0.6, COL2 * 0.6 * 0.35))

    # a: Tier1 肽长分布
    ax = axes[0]
    t1_lens = []
    for sp in novels:
        for r in novels[sp]:
            if r.get("tier") == "Tier1":
                t1_lens.append(len(r["peptide"]))
    if t1_lens:
        ax.hist(t1_lens, bins=range(7, 40, 2), color=COLORS[2], edgecolor="black", linewidth=0.3)
    ax.set_xlabel(_t("Peptide length (aa)"))
    ax.set_ylabel(_t("Count"))
    ax.set_title("a  Tier1 mirror peptides", fontsize=7, fontweight="bold", loc="left")

    # b: PSM 支持度
    ax = axes[1]
    t1_psm, t2_psm, t3_psm = [], [], []
    for sp in novels:
        for r in novels[sp]:
            tier = r.get("tier", "Tier3")
            psm = int(r.get("psm", 0))
            if tier == "Tier1":
                t1_psm.append(psm)
            elif tier == "Tier2":
                t2_psm.append(psm)
            else:
                t3_psm.append(min(psm, 50))
    data = [t for t in (t1_psm, t2_psm, t3_psm) if t]
    if data:
        bp = ax.boxplot(data, labels=["T1", "T2", "T3"], patch_artist=True)
        for patch, color in zip(bp["boxes"], COLORS[:len(data)]):
            patch.set_facecolor(color)
            patch.set_alpha(0.6)
    ax.set_ylabel(_t("PSM count (capped at 50)"))
    ax.set_title("b  PSM support by tier", fontsize=7, fontweight="bold", loc="left")
    return _save(fig, out_dir, "Fig4_mirror_validation")


# ─── 附表 ───
def tables(root, out_dir):
    try:
        import openpyxl
    except ImportError:
        return []
    wb = openpyxl.Workbook()
    # ST1: 鉴定总账
    ws1 = wb.active
    ws1.title = "ST1_identification"
    ws1.append([_t("Species"), _t("Enzyme"), _t("Tasks"), _t("Total_PSM")])
    e1 = load_tsv(os.path.join(root, "04_search_annot", "e1_report.tsv"))
    groups = defaultdict(lambda: [0, 0])
    for r in e1:
        g = groups[(r["species"], r["enzyme"])]
        g[0] += 1
        g[1] += int(r["n_psm"])
    for k in sorted(groups):
        ws1.append([k[0], k[1], groups[k][0], groups[k][1]])

    # ST2: novel 肽
    ws2 = wb.create_sheet("ST2_novel_peptides")
    ws2.append([_t("Species"), "Peptide", _t("Tier"), _t("Layer"), _t("Locus"),
                _t("Enzymes"), _t("PSMs")])
    for sp in ("Mab", "Mch"):
        for r in load_tsv(os.path.join(root, "05_novel", sp, "novel_peptides.tsv")):
            ws2.append([sp, r.get("peptide", ""), r.get("tier", ""), r.get("layer", ""),
                        r.get("locus", ""), r.get("enzymes_str", ""), r.get("psm", "")])

    # ST3: QC
    ws3 = wb.create_sheet("ST3_QC")
    ws3.append([_t("Task"), _t("Species"), _t("PSM"), _t("Peptides"),
                _t("Full_cleavage_%"), _t("MC0_%")])
    for r in load_tsv(os.path.join(root, "04a_qc", "qc_per_task.tsv")):
        ws3.append([r.get("task", ""), r.get("species", ""), r.get("n_psm", ""),
                    r.get("n_peptide", ""), r.get("pct_full_cleavage", ""), r.get("pct_mc0", "")])

    # ST4: 参数
    ws4 = wb.create_sheet("ST4_parameters")
    ws4.append([_t("Parameter"), _t("Value")])
    for k, v in [
        ("Search engine", "pFind 3.2.3"),
        ("Fixed modification", "Carbamidomethyl[C]"),
        ("Variable modifications", "Oxidation[M]; Acetyl[AnyN-term]"),
        ("Precursor tolerance", "20 ppm"),
        ("Min peptide length", "7"),
        ("Max missed cleavages", "2"),
        ("FDR", "1% (PSM & protein)"),
        ("Enzyme TL", "Trypsin KR _ C"),
        ("Enzyme LL", "LysargiNase KR _ N"),
        ("CV split", "ProteoWizard scanEvent 2/4"),
    ]:
        ws4.append([k, v])

    p = os.path.join(out_dir, "Supplementary_Tables.xlsx")
    wb.save(p)
    return [p]


# ─── Methods ───
def methods_md(root, out_dir):
    e1 = load_tsv(os.path.join(root, "04_search_annot", "e1_report.tsv"))
    e2 = load_tsv(os.path.join(root, "04_search_disc", "e2_report.tsv"))
    novels = _novels(root)
    total_psm = sum(int(r["n_psm"]) for r in e1) if e1 else 0
    total_novel = sum(len(v) for v in novels.values())
    t1 = sum(1 for v in novels.values() for r in v if r.get("tier") == "Tier1")

    # 数据集规模与物种从项目自身取，正文里不写死某个项目
    raw_rows = load_tsv(os.path.join(root, "00_input", "raw_info.tsv"))
    n_raw = len(raw_rows)
    total_gb = sum(float(r.get("size_mb") or 0) for r in raw_rows) / 1024.0
    sp_rows = load_tsv(os.path.join(root, "00_input", "species.tsv"))
    sp_names = [str(r.get("species") or r.get("species_code") or "").strip()
                for r in sp_rows]
    genome_txt = "; ".join(n for n in sp_names if n) or "the project species"

    text = f"""# Methods (auto-generated)

## Data acquisition
Raw mass spectrometry data were acquired on an Orbitrap instrument. Acquisition
parameters (instrument method, FAIMS compensation voltages, resolutions, injection
times, normalized collision energy, cycle time) are recorded per raw file in the
project raw-info table. {n_raw} raw file(s) ({total_gb:.1f} GB) were analyzed.

## Voltage splitting
Multi-compensation-voltage raw files were split by compensation voltage using
ProteoWizard msconvert with scanEvent filters derived from the acquisition method.
Peak picking used vendor-supplied centroids and MGF output carried pFind-style
titles. Split fidelity was verified by exact spectral-count reconciliation.

## Database construction
Reference genomes were obtained from NCBI RefSeq: {genome_txt}. Multi-layer
search databases were constructed: annotated proteins (L0), Prodigal-predicted
novel ORFs (L1, filtered for ≥2 novel peptides), six-frame translated ORFs
(L2, ≥10 aa, filtered), N-terminal extensions (L3), enzyme autolysis (L4),
and common contaminants (L5).

## Peptide identification
pFind (v3.2.3) was used for all searches with the following parameters:
fixed modification Carbamidomethyl[C]; variable modifications Oxidation[M]
and Acetyl[AnyN-term]; precursor tolerance 20 ppm; minimum peptide length 7;
maximum 2 missed cleavages; 1% FDR at both PSM and protein levels.
For mirror digestion samples, Trypsin+LysC (KR_C) and LysargiNase+LysN (KR_N)
enzyme specificities were applied respectively.

## Novel peptide discovery
Novel peptide candidates were identified through four-stage filtering:
(1) exclusion of peptides present in the annotated protein digestome or
identified in annotated-database searches; (2) removal of peptides that are
≥7-aa substrings of annotated proteins; (3) removal of isobaric variants
(Q↔E, GG↔N); (4) taxon-specific cleavage rules applied. Genomic locus
assignment used 200-bp neighborhood clustering.

## Mirror-digest validation
Tier1 mirror peptides were defined as novel peptides whose core sequences
(N-terminal and C-terminal K/R stripped, I/L-normalized) were identified
under both Trypsin+LysC and LysargiNase+LysN digestion conditions within
the same sample group. Tier2 required dual-enzyme support at the same
genomic locus (different peptides).

## Quantification summary
Total PSMs: {total_psm:,} | Novel peptides: {total_novel:,} |
Tier1 mirror: {t1:,}

Generated by ORFscape v{__version__} on {datetime.datetime.now().strftime('%Y-%m-%d %H:%M')}
"""
    p = os.path.join(out_dir, "Methods_auto.md")
    with open(p, "w", encoding="utf-8") as fh:
        fh.write(text)
    # 中文版 Methods（始终附一份，供内部审阅；--lang zh 时亦为主稿）
    zh = f"""# 方法（自动生成）

## 数据采集
质谱数据在 Orbitrap 仪器上采集。仪器方法、FAIMS 补偿电压、分辨率、注入时间、
归一化碰撞能量与循环时间等参数按 raw 记录在项目的 raw-info 表中。
共分析 {n_raw} 个 raw 文件（{total_gb:.1f} GB）。

## 电压拆分
多补偿电压 raw 文件使用 ProteoWizard msconvert 按采集方法推导的 scanEvent 过滤器
拆分，峰提取采用厂商质心，输出 pFind 风格标题的 MGF；拆分保真度经逐谱计数对账。

## 数据库构建
参考基因组来自 NCBI RefSeq：{genome_txt}。构建多层搜索库：注释蛋白（L0）、
Prodigal 预测新 ORF（L1）、六框翻译 ORF（L2）、N 端延伸（L3）、酶自切（L4）、
常见污染（L5）与转录证据 ORF（L6）。

## 肽段鉴定
全部搜库使用 pFind（v3.2.3）：固定修饰 Carbamidomethyl[C]；可变修饰
Oxidation[M] 与 Acetyl[AnyN-term]；前体容忍 20 ppm；最短肽长 7；最多 2 个漏切；
PSM 与蛋白水平 FDR 1%。镜像酶切样本分别按 Trypsin+LysC（KR_C）与
LysargiNase+LysN（KR_N）酶特异性搜索。

## 新肽段发现
新肽候选经四重过滤：(1) 排除注释蛋白消化集或注释库搜库已鉴定肽；(2) 剔除为注释
蛋白 ≥7 aa 连续子串的肽；(3) 剔除等重变体（Q↔E、GG↔N）；(4) 应用物种酶切规则。
基因座归并采用 200 bp 邻域聚类。

## 镜像酶切验证
Tier1 镜像肽定义为：去端 K/R 核心序列（I/L 归一）在同一样品组的 Trypsin+LysC 与
LysargiNase+LysN 双酶条件下均被鉴定的新肽。Tier2 要求同一基因座双酶支持。

## 定量摘要
总 PSM：{total_psm:,} | 新肽段：{total_novel:,} | Tier1 镜像：{t1:,}

由 ORFscape v{__version__} 生成
"""
    p_zh = os.path.join(out_dir, "Methods_auto_zh.md")
    with open(p_zh, "w", encoding="utf-8") as fh:
        fh.write(zh)
    return [p, p_zh]
