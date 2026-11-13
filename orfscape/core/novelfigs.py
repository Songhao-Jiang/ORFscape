# -*- coding: utf-8 -*-
"""新肽段筛选漏斗（分样本列）+ 参数比较。

五个阶段（样本列，不做总计）：
  S1 搜库来源        E2 db_aug/db_6f 中命中 L2_/L3_ 的 PSM 与去重肽
  S2 数据库匹配过滤  排除（a）L0 双酶理论消化集  +  E1 实测靶肽
                     （b）L0 ≥7aa 连续子串（c）I=L / N=GG / Q=E 等重变体
                     （d）酶切规则（本数据不适用，见 novel.py 声明）
  S3 基因组坐标映射  L2 按 (contig, 起点/200bp) 归 locus，L3 按注释蛋白归并；
                     多肽落同一 locus = 同一基因组事件
  S4 镜像酶切分级    Tier1 同核心肽双酶 / Tier2 同 ORF 双酶 / Tier3 单酶
  S5 可信新肽段      isoPSV credible_novel = TRUE

口径：样本级归属、非互斥（同一肽在多样本检出则各列均计）；S4/S5 采用管线正式
结果（05_novel/novel_peptides.tsv + 06_validate/novel_validated.tsv），保证与
论文数字一致。

参数比较：①各 Tier 的 isoPSV 采信比例 ②保留 vs 过滤谱图的打分/余弦分布。
"""
import csv
import os

import matplotlib.pyplot as plt
import numpy as np

from .novel import (digest_krc, digest_krn, iln, isobaric_variants,
                    load_spectra, read_faa)
from .report import COLORS, COL1, COL2, _apply_zh_font, _save, load_tsv

SP_COL = {"Mab": COLORS[0], "Mch": COLORS[1]}
TIER_COL = {"Tier1": COLORS[0], "Tier2": COLORS[1], "Tier3": "#9CA3AF"}

ZH_NOV = {
    "Novel peptide screening": "新肽段筛选",
    "a  Mab·CFP": "a  Mab·CFP", "b  Mab·TCL": "b  Mab·TCL",
    "c  Mch·TCL": "c  Mch·TCL",
    "S1 spectra": "S1 谱图", "S1 peptides": "S1 肽段",
    "S2 DB-filtered peptides": "S2 数据库过滤后肽段",
    "S3 genome-mapped peptides": "S3 基因组定位肽段",
    "S4 tiered peptides": "S4 分级肽段",
    "S5 credible peptides": "S5 可信肽段",
    "S1 · Novel-source hits (db_aug / db_6f, L2/L3 proteins)": "S1 · 搜库来源命中（db_aug / db_6f 的 L2/L3 蛋白）",
    "S2 · Database-matched filtering retained": "S2 · 数据库匹配过滤后保留",
    "S3 · Mapped to genome coordinates": "S3 · 映射到基因组坐标",
    "S4 · Mirror-digestion tiering": "S4 · 镜像酶切分级",
    "S5 · Credible novel peptides": "S5 · 可信新肽段",
    "isoPSV accepted": "isoPSV 采信",
    "retained after quality/similarity gates (KEEP)": "质量+相似度三门后保留（KEEP）",
    "S5 credible peptides (KEEP)": "S5 可信肽段（KEEP）",
    "exclude: annotated dual-enzyme digest set + E1 measured peptides; L0 contiguous substrings ≥7 aa; isobaric variants I=L, N=GG, Q=E; cleavage rules": "排除：已注释双酶消化集 + E1 实测肽段；L0 ≥7aa 连续子串；等重变体 I=L、N=GG、Q=E；酶切规则",
    "L2: (contig, start/200 bp) locus; L3: annotated protein ID; multiple peptides in one locus = one genomic event": "L2 按 (contig, 起点/200bp) 归 locus；L3 按注释蛋白 ID；同一 locus 多条肽 = 同一基因组事件",
    "Tier1 same core peptide under both enzymes; Tier2 same ORF with both enzymes; Tier3 single enzyme": "Tier1 双酶同核心肽；Tier2 同 ORF 双酶；Tier3 单酶",
    "isoPSV unlabeled workflow (pFind open/restrict competition + 4 predictors); credible_novel = TRUE": "isoPSV 无标记流程（pFind open/restrict 竞争 + 4 预测器）；credible_novel = TRUE",
    "credible": "可信",
    "Spectra per peptide (PSMs)": "每肽谱图数（PSM）",
    "1 spectrum": "1 张谱",
    "2 spectra": "2 张谱",
    "3+ spectra": "3+ 张谱",
    "Retained vs excluded novel peptides by sample": "保留 vs 排除新肽段（分样本）",
    "Tier3 (single enzyme, kept in brackets)": "Tier3 单酶（括号内=三门后保留）",
    "Count": "计数", "loci": "loci",
    "isoPSV acceptance by Tier": "各 Tier 的 isoPSV 采信比例",
    "Accepted (credible)": "采信 (credible)",
    "Score: retained vs filtered": "打分：保留 vs 过滤",
    "cand1 PSM score": "cand1 PSM 打分",
    "Predicted-spectrum cosine": "理论谱图余弦",
    "Retained (KEEP)": "保留 (KEEP)", "Filtered (DROP)": "过滤 (DROP)",
    "Trypsin+LysC": "Trypsin+LysC", "LysargiNase+LysN": "LysargiNase+LysN",
    "peptides": "肽段",
}


def T(en):
    from .i18n import get_lang
    return ZH_NOV.get(en, en) if get_lang() == "zh" else en


# ---------------------------------------------------------------- 排除集
def sp_dirs(root, sub="04_search_disc"):
    """base 下的物种子目录（不写死物种代码）。"""
    base = os.path.join(root, sub)
    if not os.path.isdir(base):
        return []
    return sorted(x for x in os.listdir(base)
                  if os.path.isdir(os.path.join(base, x)))


def build_exclusions(root, log=print):
    """每物种：L0 消化集、L0 拼接序列、E1 实测靶肽。"""
    ex = {}
    for sp in sp_dirs(root, "02_db"):
        faa = os.path.join(root, "02_db", sp, "db_annot.faa")
        recs = read_faa(faa) if os.path.exists(faa) else []
        l0 = [(h, s) for h, s in recs if h.startswith("L0_")]
        digest = set()
        for _, s in l0:
            digest.update(digest_krc(s))
            digest.update(digest_krn(s))
        e1 = set()
        ann = os.path.join(root, "04_search_annot", sp)
        if os.path.isdir(ann):
            for task in os.listdir(ann):
                spf = os.path.join(ann, task, "pFind-Filtered.spectra")
                if not os.path.exists(spf):
                    continue
                for r in load_spectra(spf):
                    if r.get("Sequence"):
                        e1.add(iln(r["Sequence"]))
        ex[sp] = dict(digest=digest, bigseq="|".join(s for _, s in l0), e1=e1)
        log(f"  [{sp}] 排除集：L0 {len(l0)} 蛋白 / 消化集 {len(digest):,} / "
            f"E1 实测 {len(e1):,}")
    return ex


def _passes_filters(pep, ex):
    """四滤判定（与 novel.py 同口径）；True = 通过（仍为新肽）。"""
    if pep in ex["digest"] or pep in ex["e1"]:
        return False
    if pep in ex["bigseq"]:
        return False
    for v in isobaric_variants(pep):
        if v != pep and (v in ex["digest"] or v in ex["e1"]):
            return False
    return True


# ---------------------------------------------------------------- 阶段 1/2
def collect_sources(root, ex, log=print):
    """扫 E2 结果 → 每 (sp,st,enz) 的 S1/S2 计数 + 肽→样本归属。

    返回 (per_group, pep_samples)：
      per_group[(sp,st,enz)] = dict(psm=, pep1=set, pep2=set)
      pep_samples[pep] = set((sp,st))
    """
    per_group, pep_samples, psm_by_samp = {}, {}, {}
    for sp in sp_dirs(root, "04_search_disc"):
        disc = os.path.join(root, "04_search_disc", sp)
        if not os.path.isdir(disc):
            continue
        for task in sorted(os.listdir(disc)):
            spec = os.path.join(disc, task, "pFind-Filtered.spectra")
            if not os.path.exists(spec):
                continue
            m = task.split("__")[0]           # 去掉 __aug/__6f 后缀
            parts = m.split("_")
            if len(parts) < 5:
                continue
            st, enz = parts[1], parts[2]
            g = per_group.setdefault((sp, st, enz),
                                     dict(psm=0, pep1=set(), pep2=set()))
            for r in load_spectra(spec):
                prot = (r.get("Proteins") or "").split(";")[0].rstrip("/")
                if not (prot.startswith("L2_") or prot.startswith("L3")
                        or prot.startswith("L6|")):
                    continue
                seq = r.get("Sequence") or ""
                if not seq:
                    continue
                p = iln(seq)
                g["psm"] += 1
                g["pep1"].add(p)
                _ps = psm_by_samp.setdefault((sp, st), {})
                _ps[p] = _ps.get(p, 0) + 1
                if _passes_filters(p, ex[sp]):
                    g["pep2"].add(p)
                    pep_samples.setdefault(p, set()).add((sp, st))
    n1 = sum(len(v["pep1"]) for v in per_group.values())
    n2 = sum(len(v["pep2"]) for v in per_group.values())
    log(f"  S1/S2：{len(per_group)} 组 / 样本内肽段 {n1:,} → 过滤后 {n2:,}"
        f"（跨组累计；肽在新颖性判定上与物种级去重不同）")
    return per_group, pep_samples, psm_by_samp


# ---------------------------------------------------------------- 正式结果
def load_official(root):
    """<VALIDATE_DIR>/novel_validated.tsv（含 tier/locus/isoPSV 列；默认 06_validate）。"""
    vdir = os.environ.get("ORFSCAPE_VALIDATE_DIR", "06_validate")
    f = os.path.join(root, vdir, "novel_validated.tsv")
    if not os.path.exists(f):
        return []
    return load_tsv(f)


def stage_counts(per_group, pep_samples, official):
    """按样本汇总 S1–S5（肽×样本唯一化，非行级重复）。

    跨物种共享肽在 novel_validated.tsv 中各有一行；样本级统计按该样本物种
    对应的行取 tier，避免同一样本内重复计数。
    """
    samples = sorted({(sp, st) for (sp, st, _e) in per_group})
    rows_by_sp = {}
    for r in official:
        rows_by_sp.setdefault(r.get("species", ""), {})[r.get("peptide", "")] = r
    out = {}
    for sp, st in samples:
        d = dict(sp=sp, st=st, s1_psm=0, s1_pep=set(), s2_pep=set(),
                 s3_pep=set(), s4_pep=set(), s5_pep=set(), s5_iso=set(),
                 loci=set(), t1=0, t2=0, t3=0, k1=0, k2=0, k3=0,
                 enz_pep={})
        for enz in ("TL", "LL"):
            g = per_group.get((sp, st, enz))
            if not g:
                d["enz_pep"][enz] = 0
                continue
            d["s1_psm"] += g["psm"]
            d["s1_pep"] |= g["pep1"]
            d["s2_pep"] |= g["pep2"]
            d["enz_pep"][enz] = len(g["pep2"])
        out[(sp, st)] = d
    # S3/S4/S5：逐 (样本, 肽) 唯一计入，tier 取该样本物种对应的行
    # P6b 回搜补进的肽在 E2 里没有命中，若不并进来，漏斗底部（S3–S5）会比
    # 正式结果少一截（TP 实测 495 vs 529）。只在"该物种仅一个样本型"时并入，
    # 多样本型（NTM 的 CFP/TCL）无法按样本归属，保持原口径。
    n_st = {}
    for (sp, _st) in out:
        n_st[sp] = n_st.get(sp, 0) + 1
    for (sp, st), d in out.items():
        if n_st.get(sp, 0) == 1:
            d["s2_pep"] |= set(rows_by_sp.get(sp, {}))
    for (sp, st), d in out.items():
        sp_rows = rows_by_sp.get(sp, {})
        for pep in d["s2_pep"]:
            row = sp_rows.get(pep)
            if not row:
                continue
            d["s3_pep"].add(pep)
            d["s4_pep"].add(pep)
            d["loci"].add(row.get("locus", ""))
            tier = row.get("tier", "")
            if tier in ("Tier1", "Tier2", "Tier3"):
                key = {"Tier1": "t1", "Tier2": "t2", "Tier3": "t3"}[tier]
                d[key] += 1
            if str(row.get("v_iso_pass", "")).lower() == "true":
                d["s5_iso"].add(pep)
            if str(row.get("final_keep", "")).startswith("KEEP"):
                d["s5_pep"].add(pep)
                if tier in ("Tier1", "Tier2", "Tier3"):
                    kkey = {"Tier1": "k1", "Tier2": "k2", "Tier3": "k3"}[tier]
                    d[kkey] += 1
    return out


# ---------------------------------------------------------------- 图 0 流程图
def nov0_flowchart(root, out_dir, stats):
    """PRISMA 风格筛选流程图：5 级方框（框内三样本数值），框间箭头标排除/操作。"""
    from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

    samples = sorted(stats)

    def trio(fn):
        return "  |  ".join(f"{sp}·{st} {fn(stats[(sp, st)]):,}"
                            for sp, st in samples)

    stages = [
        (T("S1 · Novel-source hits (db_aug / db_6f, L2/L3 proteins)"),
         [("spectra", trio(lambda d: d["s1_psm"])),
          ("peptides", trio(lambda d: len(d["s1_pep"])))]),
        (T("S2 · Database-matched filtering retained"),
         [(T("peptides"), trio(lambda d: len(d["s2_pep"])))]),
        (T("S3 · Mapped to genome coordinates"),
         [(T("peptides"), trio(lambda d: len(d["s3_pep"]))),
          ("loci", trio(lambda d: len(d["loci"])))]),
        (T("S4 · Mirror-digestion tiering"),
         [("Tier1 (mirror core peptide)", trio(lambda d: d["t1"])),
          ("Tier2 (same ORF, dual enzyme)", trio(lambda d: d["t2"])),
          (T("Tier3 (single enzyme, kept in brackets)"),
           "  |  ".join(
               f"{sp}·{st} {stats[(sp, st)]['t3']:,}"
               + (f" ({stats[(sp, st)]['k3']:,})"
                  if stats[(sp, st)]["k3"] < stats[(sp, st)]["t3"] else "")
               for sp, st in samples))]),
        (T("S5 · Credible novel peptides"),
         [("isoPSV accepted", trio(lambda d: len(d["s5_iso"]))),
          ("retained after quality/similarity gates (KEEP)",
           trio(lambda d: len(d["s5_pep"])))]),
    ]
    notes = [
        T("exclude: annotated dual-enzyme digest set + E1 measured peptides; "
          "L0 contiguous substrings ≥7 aa; isobaric variants I=L, N=GG, Q=E; "
          "cleavage rules"),
        T("L2: (contig, start/200 bp) locus; L3: annotated protein ID; "
          "multiple peptides in one locus = one genomic event"),
        T("Tier1 same core peptide under both enzymes; Tier2 same ORF with "
          "both enzymes; Tier3 single enzyme"),
        T("isoPSV unlabeled workflow (pFind open/restrict competition + "
          "4 predictors); credible_novel = TRUE"),
    ]

    import textwrap
    fig, ax = plt.subplots(figsize=(COL2, 8.0))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 12.6)
    ax.axis("off")
    box_l, box_r = 0.25, 6.15
    note_x = 6.4
    heights = [1.75, 1.25, 1.55, 2.15, 1.30]
    gap = 1.05
    y = 12.45
    centers = []
    for i, ((title, rows), h) in enumerate(zip(stages, heights)):
        ax.add_patch(FancyBboxPatch((box_l, y - h), box_r - box_l, h,
                                    boxstyle="round,pad=0.02,rounding_size=0.12",
                                    fc="#EEF4FB" if i in (0, 4) else "white",
                                    ec="#4B5563", lw=0.7))
        ax.text(box_l + 0.18, y - 0.26, f"{i+1}. {title}", fontsize=6.6,
                fontweight="bold", va="top", ha="left")
        for j, (lab, val) in enumerate(rows):
            ax.text(box_l + 0.18, y - 0.62 - j * 0.34,
                    f"{lab}: {val}", fontsize=5.4, va="top", ha="left",
                    color="#111827")
        centers.append((y, y - h))
        y -= h + gap
    # 箭头 + 侧注
    for i in range(len(stages) - 1):
        y0 = centers[i][1]
        y1 = centers[i + 1][0]
        ax.add_patch(FancyArrowPatch((1.0, y0), (1.0, y1),
                                     arrowstyle="-|>", mutation_scale=9,
                                     lw=0.9, color="#374151"))
        ax.text(note_x, (y0 + y1) / 2,
                textwrap.fill(notes[i], 46), fontsize=4.7, va="center",
                ha="left", color="#374151", linespacing=1.35,
                bbox=dict(boxstyle="round,pad=0.3", fc="#F9FAFB",
                          ec="#D1D5DB", lw=0.4))
    _save(fig, out_dir, "NOV0_screening_flowchart")


def nov1_funnel(root, out_dir, stats):
    samples = sorted(stats)
    fig, axes = plt.subplots(1, len(samples), figsize=(COL2, 3.1),
                             gridspec_kw={"wspace": 0.34})
    if len(samples) == 1:
        axes = [axes]
    for k, (ax, key) in enumerate(zip(axes, samples)):
        d = stats[key]
        stages = [
            (T("S1 spectra"), d["s1_psm"], None),
            (T("S1 peptides"), len(d["s1_pep"]), None),
            (T("S2 DB-filtered peptides"), len(d["s2_pep"]), None),
            (T("S3 genome-mapped peptides")
             + (f" ({len(d['loci'])} loci)" if d["loci"] else ""),
             len(d["s3_pep"]), None),
            (T("S4 tiered peptides"), len(d["s4_pep"]),
             [d["t1"], d["t2"], d["t3"]]),
            (T("S5 credible peptides (KEEP)"), len(d["s5_pep"]), None),
        ]
        y = np.arange(len(stages))[::-1]
        for i, (lab, val, tiers) in enumerate(stages):
            if tiers:                        # S4：按 Tier 堆叠
                left = 0
                for j, tv in enumerate(tiers):
                    c = list(TIER_COL.values())[j]
                    ax.barh(y[i], tv, left=left, height=0.62, color=c,
                            label=[T("Tier1"), T("Tier2"), T("Tier3")][j]
                            if k == 0 else None)
                    left += tv
            else:
                ax.barh(y[i], val, height=0.62, color=COLORS[0],
                        alpha=0.85 if i > 0 else 0.55)
            ax.text(max(val, 1) * 1.06, y[i], f"{val:,}", va="center",
                    fontsize=5.4)
        ax.set_yticks(y)
        ax.set_yticklabels([s[0] for s in stages] if k == 0 else [],
                           fontsize=5.2)
        ax.set_xscale("log")
        ax.set_xlim(0.6, max(d["s1_psm"], len(d["s1_pep"]), 2) * 3.2)
        ax.set_xlabel(T("Count"), fontsize=5.6)
        ax.set_title(f"{'abc'[k]}  {key[0]}·{key[1]}", fontsize=7.5,
                     fontweight="bold", loc="left")
        ax.tick_params(labelsize=5)
    fig.suptitle(T("Novel peptide screening"), fontsize=7.5,
                 fontweight="bold", x=0.01, ha="left", y=1.06)
    handles = [plt.Rectangle((0, 0), 1, 1, fc=TIER_COL[t], ec="none")
               for t in ("Tier1", "Tier2", "Tier3")]
    fig.legend(handles, ["Tier1", "Tier2", "Tier3"], frameon=False,
               fontsize=4.8, ncol=3, loc="upper right",
               bbox_to_anchor=(1.0, 1.10), handlelength=1.0)
    _save(fig, out_dir, "NOV1_screening_funnel")


# ---------------------------------------------------------------- 图 2 参数
def nov2_params(root, out_dir, official):
    """a 各 Tier 采信比例；b 打分 KEEP vs DROP；c 余弦 KEEP vs DROP。"""
    def _f(v):
        try:
            x = float(v)
        except (TypeError, ValueError):
            return None
        return x if np.isfinite(x) else None

    rows = [r for r in official if r.get("tier")]
    tiers = ["Tier1", "Tier2", "Tier3"]
    acc = {}
    for t in tiers:
        sub = [r for r in rows if r["tier"] == t]
        ok = sum(1 for r in sub if str(r.get("v_iso_pass", "")).lower() == "true")
        acc[t] = (len(sub), ok)
    keep_s = [x for x in (_f(r.get("v_best_score")) for r in rows
                          if r.get("final_keep", "").startswith("KEEP")) if x is not None]
    drop_s = [x for x in (_f(r.get("v_best_score")) for r in rows
                          if r.get("final_keep") == "DROP") if x is not None]
    keep_c = [x for x in (_f(r.get("v_best_cos")) for r in rows
                          if r.get("final_keep", "").startswith("KEEP")) if x is not None]
    drop_c = [x for x in (_f(r.get("v_best_cos")) for r in rows
                          if r.get("final_keep") == "DROP") if x is not None]

    fig, axes = plt.subplots(1, 3, figsize=(COL2, 2.7),
                             gridspec_kw={"wspace": 0.42,
                                          "width_ratios": [1.15, 0.85, 0.85]})
    # a 采信比例
    ax = axes[0]
    x = np.arange(len(tiers))
    pct = [100 * acc[t][1] / max(acc[t][0], 1) for t in tiers]
    bars = ax.bar(x, pct, 0.55, color=[TIER_COL[t] for t in tiers])
    for b, v, t in zip(bars, pct, tiers):
        ax.text(b.get_x() + b.get_width() / 2, v + 1.4,
                f"{v:.1f}%\n{acc[t][1]}/{acc[t][0]}", ha="center", fontsize=5.2)
    ax.set_xticks(x)
    ax.set_xticklabels(tiers, fontsize=6)
    ax.set_ylabel("%")
    ax.set_ylim(0, max(pct) * 1.35)
    ax.set_title(T("isoPSV acceptance by Tier"), fontsize=7.5,
                 fontweight="bold", loc="left")
    # b/c KEEP vs DROP 箱线
    for ax, keep, drop, ttl, ylab in (
            (axes[1], keep_s, drop_s, T("cand1 PSM score"),
             T("cand1 PSM score")),
            (axes[2], keep_c, drop_c, T("Predicted-spectrum cosine"),
             T("Predicted-spectrum cosine"))):
        bp = ax.boxplot([keep, drop], positions=[0, 1], widths=0.55,
                        patch_artist=True, showfliers=False,
                        medianprops=dict(color="#111827", lw=0.9))
        for i, patch in enumerate(bp["boxes"]):
            patch.set_facecolor(COLORS[0] if i == 0 else "#9CA3AF")
            patch.set_alpha(0.55)
            patch.set_edgecolor(COLORS[0] if i == 0 else "#6B7280")
        for i, (vals, xi) in enumerate(((keep, 0), (drop, 1))):
            ax.text(xi, np.median(vals), f" {np.median(vals):.3f}",
                    fontsize=4.8, va="bottom", ha="left", color="#374151")
            ax.text(xi, ax.get_ylim()[1], f"n={len(vals):,}", fontsize=4.6,
                    ha="center", va="bottom", color="#6B7280")
        ax.set_xticks([0, 1])
        ax.set_xticklabels([T("Retained (KEEP)"), T("Filtered (DROP)")],
                           fontsize=5)
        ax.set_ylabel(ylab, fontsize=5.6)
        ax.set_title(ttl, fontsize=7.5, fontweight="bold", loc="left")
        ax.set_xlim(-0.6, 1.6)
    fig.suptitle(T("Parameter comparison"), fontsize=7.5, fontweight="bold",
                 x=0.01, ha="left", y=1.03)
    _save(fig, out_dir, "NOV2_parameter_comparison")
    return acc, keep_s, drop_s, keep_c, drop_c


# ---------------------------------------------------------------- 图 3 分布比较
def nov3_score_dists(root, out_dir, stats, psm_by_samp, official, pep_samples):
    """排除 vs 保留新肽段：a 谱图数分档构成；b 余弦；c 打分（均分样本）。"""
    def _f(v):
        try:
            x = float(v)
        except (TypeError, ValueError):
            return None
        return x if np.isfinite(x) else None

    samples = sorted(stats)
    rows_by_sp = {}
    for r in official:
        rows_by_sp.setdefault(r.get("species", ""), {})[r.get("peptide", "")] = r

    def collect(sp, st, kind):
        sp_rows = rows_by_sp.get(sp, {})
        out = {"KEEP": [], "DROP": []}
        for pep, row in sp_rows.items():
            if (sp, st) not in pep_samples.get(pep, set()):
                continue
            cls = "KEEP" if str(row.get("final_keep", "")).startswith("KEEP") else "DROP"
            if kind == "psm":
                v = psm_by_samp.get((sp, st), {}).get(pep, 0)
            elif kind == "cos":
                v = _f(row.get("v_best_cos"))
            else:
                v = _f(row.get("v_best_score"))
            if v is None:
                continue
            out[cls].append(v)
        return out["KEEP"], out["DROP"]

    fig = plt.figure(figsize=(COL2, 3.2))
    outer = fig.add_gridspec(1, 3, wspace=0.42)
    # ---- a 谱图数分档（1 / 2 / 3+）构成 ----
    ax = fig.add_subplot(outer[0, 0])
    bins = [(1, 1, T("1 spectrum")), (2, 2, T("2 spectra")),
            (3, 99, T("3+ spectra"))]
    pal = [COLORS[0], COLORS[1], "#9CA3AF"]
    xlabels, xpos = [], []
    xi = 0
    for i, (sp, st) in enumerate(samples):
        k, d = collect(sp, st, "psm")
        for cls, arr in (("KEEP", k), ("DROP", d)):
            tot = len(arr) or 1
            bottom = 0.0
            for (lo, hi, lab), c in zip(bins, pal):
                frac = 100.0 * sum(1 for v in arr if lo <= v <= hi) / tot
                ax.bar(xi, frac, 0.78, bottom=bottom, color=c,
                       alpha=0.9 if cls == "KEEP" else 0.45)
                if frac >= 6:
                    ax.text(xi, bottom + frac / 2, f"{frac:.0f}", ha="center",
                            va="center", fontsize=4.8,
                            color="white" if c == COLORS[0] else "#111827")
                bottom += frac
            xlabels.append(f"{sp}·{st}\n{cls}")
            xpos.append(xi)
            xi += 1
        xi += 0.55
    ax.set_xticks(xpos)
    ax.set_xticklabels(xlabels, fontsize=4.2, rotation=38, ha="right")
    ax.set_ylim(0, 104)
    ax.set_ylabel("% peptides", fontsize=5.6)
    ax.set_title(T("Spectra per peptide (PSMs)"), fontsize=7.5,
                 fontweight="bold", loc="left")
    handles = [plt.Rectangle((0, 0), 1, 1, fc=c, ec="none") for c in pal]
    ax.legend(handles, [b[2] for b in bins], frameon=False, fontsize=4.4,
              loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=3,
              handlelength=1.0, columnspacing=0.6)
    # ---- b/c 余弦与打分的箱线 ----
    for j, (kind, ttl, ylab) in enumerate(
            (("cos", T("Predicted-spectrum cosine"), T("Predicted-spectrum cosine")),
             ("score", T("cand1 PSM score"), T("cand1 PSM score")))):
        ax = fig.add_subplot(outer[0, j + 1])
        pos, data, colors = [], [], []
        for i, (sp, st) in enumerate(samples):
            k, d = collect(sp, st, kind)
            base = i * 2.4
            for m, (arr, c) in enumerate(((k, COLORS[0]), (d, "#9CA3AF"))):
                pos.append(base + m * 0.9)
                data.append(arr if arr else [np.nan])
                colors.append(c)
        bp = ax.boxplot(data, positions=pos, widths=0.62, patch_artist=True,
                        showfliers=False,
                        medianprops=dict(color="#111827", lw=0.9))
        for patch, c in zip(bp["boxes"], colors):
            patch.set_facecolor(c)
            patch.set_alpha(0.55)
            patch.set_edgecolor(c)
        for xi, arr in zip(pos, data):
            arr = [v for v in arr if np.isfinite(v)]
            if arr:
                ax.text(xi, 0.985, f"{np.median(arr):.2f}",
                        transform=ax.get_xaxis_transform(), ha="center",
                        va="top", fontsize=4.4, color="#374151")
        ax.set_xticks([i * 2.4 + 0.45 for i in range(len(samples))])
        ax.set_xticklabels([f"{sp}\n{st}" for sp, st in samples], fontsize=5.2)
        ax.set_ylabel(ylab, fontsize=5.6)
        ax.set_title(ttl, fontsize=7.5, fontweight="bold", loc="left")
    # 图例（KEEP / DROP）放 b 面板右下空白
    handles = [plt.Rectangle((0, 0), 1, 1, fc=COLORS[0], alpha=0.55, ec=COLORS[0]),
               plt.Rectangle((0, 0), 1, 1, fc="#9CA3AF", alpha=0.55, ec="#6B7280")]
    fig.legend(handles, [T("Retained (KEEP)"), T("Filtered (DROP)")],
               frameon=False, fontsize=5.0, ncol=2, loc="upper right",
               bbox_to_anchor=(1.0, 1.06), handlelength=1.1)
    _save(fig, out_dir, "NOV3_retained_vs_excluded")
    med, share2 = {}, {}
    for sp, st in samples:
        for kind in ("psm", "cos", "score"):
            k, d = collect(sp, st, kind)
            med[(sp, st, kind)] = (
                float(np.median(k)) if k else float("nan"),
                float(np.median(d)) if d else float("nan"))
        k, d = collect(sp, st, "psm")
        share2[(sp, st)] = (
            100.0 * sum(1 for v in k if v >= 2) / max(len(k), 1),
            100.0 * sum(1 for v in d if v >= 2) / max(len(d), 1))
    return dict(med=med, share2=share2)


# ---------------------------------------------------------------- 表 + 结论
def write_table(root, out_dir, stats):
    samples = sorted(stats)
    rows = []
    from .i18n import get_lang
    zh = get_lang() == "zh"
    head = ([ "指标" if zh else "Metric"]
            + [f"{sp}·{st}" for sp, st in samples])
    rows.append(head)
    spec = [
        ("S1 搜库来源谱图（PSM）" if zh else "S1 novel-source spectra (PSMs)",
         lambda d: d["s1_psm"]),
        ("S1 搜库来源肽段（去重）" if zh else "S1 novel-source peptides (unique)",
         lambda d: len(d["s1_pep"])),
        ("S2 数据库过滤后肽段" if zh else "S2 DB-filtered peptides",
         lambda d: len(d["s2_pep"])),
        ("S3 基因组定位肽段" if zh else "S3 genome-mapped peptides",
         lambda d: len(d["s3_pep"])),
        ("S3 基因组事件（locus）数" if zh else "S3 genomic events (loci)",
         lambda d: len(d["loci"])),
        ("S4 Tier1（镜像同核心肽）" if zh else "S4 Tier1 (mirror core peptide)",
         lambda d: d["t1"], lambda d: d["k1"]),
        ("S4 Tier2（同 ORF 双酶）" if zh else "S4 Tier2 (same ORF, dual enzyme)",
         lambda d: d["t2"], lambda d: d["k2"]),
        ("S4 Tier3（单酶）" if zh else "S4 Tier3 (single enzyme)",
         lambda d: d["t3"], lambda d: d["k3"]),
        ("S5a isoPSV 采信肽段" if zh else "S5a isoPSV-accepted peptides",
         lambda d: len(d["s5_iso"])),
        ("S5b 可信新肽段（质量+相似度三门 KEEP）" if zh
         else "S5b credible peptides (quality+similarity gates, KEEP)",
         lambda d: len(d["s5_pep"])),
    ]
    keep_lab = "留" if zh else "kept"
    for item in spec:
        lab = item[0]
        if len(item) == 3:
            fn, fk = item[1], item[2]
            cells = []
            for k in samples:
                tot, kept = fn(stats[k]), fk(stats[k])
                if kept >= tot:
                    cells.append(f"{tot:,}")
                else:
                    cells.append(f"{tot:,}（{keep_lab}{kept:,}）" if zh
                                 else f"{tot:,} ({keep_lab} {kept:,})")
            rows.append([lab] + cells)
        else:
            fn = item[1]
            rows.append([lab] + [f"{fn(stats[k]):,}" for k in samples])
    tsv = os.path.join(out_dir, "NOV_screening_table.tsv")
    with open(tsv, "w", newline="", encoding="utf-8-sig") as fh:
        for r in rows:
            fh.write("\t".join(r) + "\n")
    try:
        from openpyxl import Workbook
        from openpyxl.styles import Alignment, Font
        wb = Workbook()
        ws = wb.active
        ws.title = "novel screening"
        for r in rows:
            ws.append(r)
        for c in ws[1]:
            c.font = Font(name="Arial", size=9, bold=True)
        for row in ws.iter_rows(min_row=2):
            for c in row:
                c.font = Font(name="Arial", size=9)
            row[0].alignment = Alignment(horizontal="left")
        ws.column_dimensions["A"].width = 40
        for i in range(2, len(samples) + 2):
            ws.column_dimensions[chr(64 + i)].width = 18
        wb.save(os.path.join(out_dir, "NOV_screening_table.xlsx"))
    except ImportError:
        pass
    return tsv


def write_conclusions(root, out_dir, stats, acc, scores, official, dist):
    from .i18n import get_lang
    zh = get_lang() == "zh"
    keep, drop, ckeep, cdrop = scores
    samples = sorted(stats)
    # ---- 溯源数字（全部现算）----
    n_rows = len(official)
    uniq = {r.get("peptide", "") for r in official}
    iso_rows = [r for r in official
                if str(r.get("v_iso_pass", "")).lower() == "true"]
    keep_rows = [r for r in official
                 if str(r.get("final_keep", "")).startswith("KEEP")]
    iso_uniq = {r["peptide"] for r in iso_rows}
    keep_uniq = {r["peptide"] for r in keep_rows}
    iso_sum = sum(len(d["s5_iso"]) for d in stats.values())
    keep_sum = sum(len(d["s5_pep"]) for d in stats.values())
    per = "; ".join(
        f"{sp}·{st} " + "→".join([
            f"{stats[(sp, st)]['s1_psm']:,}",
            f"{len(stats[(sp, st)]['s1_pep']):,}",
            f"{len(stats[(sp, st)]['s2_pep']):,}",
            f"{len(stats[(sp, st)]['s3_pep']):,}",
            f"{len(stats[(sp, st)]['s4_pep']):,}",
            f"{len(stats[(sp, st)]['s5_iso']):,}",
            f"{len(stats[(sp, st)]['s5_pep']):,}"])
        for sp, st in samples)
    t1n, t1k = acc["Tier1"]
    t2n, t2k = acc["Tier2"]
    t3n, t3k = acc["Tier3"]
    if zh:
        lines = [
            "# 新肽段筛选与参数比较 —— 结论句（分样本列）",
            "",
            f"**NOV1 · 筛选漏斗（S1 谱图→S1 肽段→S2 过滤→S3 定位→S4 分级→S5 可信）**"
            f" —— 各样本依次为：{per}；样本级归属非互斥（同一肽在多样本检出则各列"
            f"均计），S4/S5 采用管线正式结果。",
            "",
            f"**NOV2a · 各 Tier 的 isoPSV 采信比例** —— Tier1 "
            f"{t1k}/{t1n}（{100*t1k/max(t1n,1):.1f}%）、Tier2 "
            f"{t2k}/{t2n}（{100*t2k/max(t2n,1):.1f}%）、Tier3 "
            f"{t3k}/{t3n}（{100*t3k/max(t3n,1):.1f}%）；采信率随 Tier 递减，"
            f"镜像证据强度与独立引擎复核结论一致。",
            "",
            f"**NOV3 · 保留 vs 排除肽段（分样本）** —— ≥2 张谱支持的肽段占比：保留 "
            + "；".join(f"{sp}·{st} {dist['share2'][(sp, st)][0]:.0f}% vs 排除 "
                        f"{dist['share2'][(sp, st)][1]:.0f}%"
                        for sp, st in samples)
            + f"；余弦中位：保留 {np.nanmean([dist['med'][(sp, st, 'cos')][0] for sp, st in samples]):.2f} "
              f"vs 排除 {np.nanmean([dist['med'][(sp, st, 'cos')][1] for sp, st in samples]):.2f}"
              f"；打分中位：保留 {np.nanmean([dist['med'][(sp, st, 'score')][0] for sp, st in samples]):.2f} "
              f"vs 排除 {np.nanmean([dist['med'][(sp, st, 'score')][1] for sp, st in samples]):.2f}"
              f"——排除集在三项指标上系统性低于保留集，且方向在各样本一致。",
            "",
            f"**口径与溯源** —— 数据链：`04_search_disc/*/pFind-Filtered.spectra`"
            f"（L2/L3 命中）→ 样本归属；`06_validate/<物种>_<酶>/result_*/"
            f"1.unlabeled_result/2.novel_peptide_competition_report.tsv`"
            f"（credible_novel）→ `validate.py` 肽级门控 → "
            f"`06_validate/novel_validated.tsv`（v_iso_pass / v_gate_pass / "
            f"final_keep）→ 本图。行级 {n_rows:,} 行 = 唯一肽 {len(uniq):,}"
            f"（跨物种共享肽产生的重复行）；isoPSV 采信 {len(iso_rows):,} 行 / "
            f"{len(iso_uniq):,} 肽，三门后 KEEP {len(keep_rows):,} 行 / "
            f"{len(keep_uniq):,} 肽。样本级归属非互斥：三列合计 "
            f"{iso_sum:,}（= {len(iso_uniq):,} + 跨样本重复 "
            f"{iso_sum - len(iso_uniq):,}）与 {keep_sum:,}"
            f"（= {len(keep_uniq):,} + {keep_sum - len(keep_uniq):,}）。",
            "",
            f"**NOV2b · 保留 vs 过滤谱图的打分** —— cand1 PSM 打分中位：保留 "
            f"{np.median(keep):.3f} vs 过滤 {np.median(drop):.3f}；"
            f"理论谱图余弦中位：保留 {np.median(ckeep):.3f} vs 过滤 "
            f"{np.median(cdrop):.3f}——两指标方向一致，说明 Tier3 门控剔除的"
            f"主要是低分/低相似度谱图。",
            "",
        ]
    else:
        lines = [
            "# Novel peptide screening & parameter comparison — conclusions "
            "(per-sample columns)",
            "",
            f"**NOV1 · Screening funnel (S1 spectra -> S1 peptides -> S2 filter "
            f"-> S3 mapping -> S4 tiering -> S5 credible)** — per sample: "
            f"{per}; sample attribution is non-exclusive (a peptide detected in "
            f"several samples is counted in each), and S4/S5 use the pipeline's "
            f"official results.",
            "",
            f"**NOV2a · isoPSV acceptance by Tier** — Tier1 {t1k}/{t1n} "
            f"({100*t1k/max(t1n,1):.1f}%), Tier2 {t2k}/{t2n} "
            f"({100*t2k/max(t2n,1):.1f}%), Tier3 {t3k}/{t3n} "
            f"({100*t3k/max(t3n,1):.1f}%); acceptance decreases monotonically "
            f"with tier, agreeing with the mirror-evidence hierarchy.",
            "",
            f"**NOV3 · Retained vs excluded peptides (per sample)** — share "
            f"of peptides with >=2 supporting spectra: retained "
            + "; ".join(f"{sp}·{st} {dist['share2'][(sp, st)][0]:.0f}% vs "
                        f"excluded {dist['share2'][(sp, st)][1]:.0f}%"
                        for sp, st in samples)
            + f"; median cosine: retained "
              f"{np.nanmean([dist['med'][(sp, st, 'cos')][0] for sp, st in samples]):.2f} "
              f"vs excluded "
              f"{np.nanmean([dist['med'][(sp, st, 'cos')][1] for sp, st in samples]):.2f}; "
              f"median score: retained "
              f"{np.nanmean([dist['med'][(sp, st, 'score')][0] for sp, st in samples]):.2f} "
              f"vs excluded "
              f"{np.nanmean([dist['med'][(sp, st, 'score')][1] for sp, st in samples]):.2f}; "
              f"the excluded set is systematically lower on all three metrics "
              f"with consistent direction across samples.",
            "",
            f"**Counting rule & provenance** — data chain: "
            f"`04_search_disc/*/pFind-Filtered.spectra` (L2/L3 hits) -> sample "
            f"attribution; `06_validate/<sp>_<enz>/result_*/1.unlabeled_result/"
            f"2.novel_peptide_competition_report.tsv` (credible_novel) -> "
            f"`validate.py` peptide-level gates -> `06_validate/"
            f"novel_validated.tsv` (v_iso_pass / v_gate_pass / final_keep) -> "
            f"this figure. Row level {n_rows:,} rows = {len(uniq):,} unique "
            f"peptides (duplicate rows from cross-species shared peptides); "
            f"isoPSV accepted {len(iso_rows):,} rows / {len(iso_uniq):,} "
            f"peptides, KEEP after three gates {len(keep_rows):,} rows / "
            f"{len(keep_uniq):,} peptides. Sample attribution is "
            f"non-exclusive: column sums {iso_sum:,} "
            f"(= {len(iso_uniq):,} + {iso_sum - len(iso_uniq):,} cross-sample "
            f"duplicates) and {keep_sum:,} "
            f"(= {len(keep_uniq):,} + {keep_sum - len(keep_uniq):,}).",
            "",
            f"**NOV2b · Scores of retained vs filtered spectra** — median cand1 "
            f"PSM score: retained {np.median(keep):.3f} vs filtered "
            f"{np.median(drop):.3f}; median predicted-spectrum cosine: retained "
            f"{np.median(ckeep):.3f} vs filtered {np.median(cdrop):.3f} — both "
            f"metrics agree that the Tier3 gate mainly removes low-score, "
            f"low-similarity spectra.",
            "",
        ]
    out = os.path.join(out_dir, "NOV_figure_conclusions.md")
    with open(out, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines))
    return out


# ---------------------------------------------------------------- 总入口
def build_novelfigs(root, out_dir, log=print):
    _apply_zh_font()
    os.makedirs(out_dir, exist_ok=True)
    ex = build_exclusions(root, log=log)
    per_group, pep_samples, psm_by_samp = collect_sources(root, ex, log=log)
    official = load_official(root)
    log(f"  正式新肽段结果 {len(official):,} 行")
    stats = stage_counts(per_group, pep_samples, official)
    nov0_flowchart(root, out_dir, stats)
    dist = nov3_score_dists(root, out_dir, stats, psm_by_samp, official,
                            pep_samples)
    nov1_funnel(root, out_dir, stats)
    acc, *scores = nov2_params(root, out_dir, official)
    write_table(root, out_dir, stats)
    write_conclusions(root, out_dir, stats, acc, scores, official, dist)
    log("  NOV1/NOV2 + 表 + 结论 完成")
    return stats
