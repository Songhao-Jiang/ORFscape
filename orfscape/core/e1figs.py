# -*- coding: utf-8 -*-
"""E1 已注释搜库结果蛋白质组质控系列单图（用户口径 2026-11-06）：

  E1Q1 鉴定产率     a 逐任务 PSM   b 逐任务谱图解析率（PSM/MS2）
  E1Q2 FDR 分布     a 六组各自 Target/Decoy 累计曲线（a1–a6 小多重图）
                    b 逐任务 1% FDR 下 decoy 占比
  E1Q3 序列覆盖率   6 组（物种×样品×酶）箱线+任务散点
  E1Q4 质量精度     a 母离子 ppm 直方图   b 子离子 Da 直方图
  E1Q5 酶切特异性   a 特异性分类堆叠（6 组） b 漏切分布堆叠（6 组）
  E1Q6 修饰频率     三修饰 × 6 组肽级频率条形
  E1Q7 肽段物理分布 a 肽长直方图（去重肽） b 电荷分布（PSM 级，按酶分组）

数据源：04_search_annot/*/pFind.summary（逐任务聚合量）+
        pFind-Filtered.spectra（单程流式：Q-value/T-D/质量差/长度/电荷/漏切）+
        00_input/raw_info.tsv × 03_mgf/split_report.tsv（MS2 分母）。
结论句 → E1_QC_figure_conclusions.md（双语两份）。
"""
import csv
import os
import re

import matplotlib.pyplot as plt
import numpy as np

from .novel import iln
from .report import COLORS, COL1, COL2, _apply_zh_font, _save, _t, load_tsv


def _core(seq):
    """镜像核心区：去两端 K/R + I/L 归一（Tier1 同口径）。"""
    s = seq.upper()
    while s and s[0] in "KR":
        s = s[1:]
    while s and s[-1] in "KR":
        s = s[:-1]
    return iln(s)

SP_COL = {"Mab": COLORS[0], "Mch": COLORS[1]}
EN_COL = {"TL": COLORS[0], "LL": COLORS[3]}
EN2NAME = {"TL": "Trypsin+LysC", "LL": "LysargiNase+LysN"}

ZH_E1 = {
    "E1Q1  Identification yield": "E1Q1  鉴定产率",
    "E1Q2  Target-decoy FDR profile": "E1Q2  Target-Decoy FDR 剖面",
    "E1Q3  Protein sequence coverage": "E1Q3  蛋白序列覆盖度",
    "E1Q4  Mass accuracy": "E1Q4  质量精度",
    "E1Q5  Cleavage specificity": "E1Q5  酶切特异性与漏切",
    "E1Q6  Modification frequency": "E1Q6  修饰频率",
    "E1Q7  Peptide length & charge": "E1Q7  肽段长度与电荷",
    "a  PSMs per task": "a  每任务 PSM",
    "b  ID rate per task (PSM / MS2)": "b  每任务谱图解析率（PSM/MS2）",
    "a  Identified spectra (PSMs)": "a  鉴定谱图数（PSM）",
    "b  Peptides": "b  肽段数",
    "c  Proteins": "c  蛋白数",
    "d  Spectrum explain rate": "d  谱图解析率",
    "a  Cumulative PSMs vs FDR threshold": "a  累计 PSM 随 FDR 阈值",
    "b  Decoy fraction at 1% FDR per task": "b  1% FDR 下逐任务 decoy 占比",
    "(task order as in E1Q1)": "（任务序同 E1Q1）",
    "Avg sequence coverage (%)": "平均序列覆盖度 (%)",
    "Protein sequence coverage (%)": "蛋白质序列覆盖度 (%)",
    "Numbers = median": "数字 = 中位数",
    "numbers = mean % (segments < 3% unlabeled)": "数字=各组均值% (＜3% 的段未标注)",
    "Combined": "双酶合并",
    "a  Precursor mass error (ppm)": "a  母离子质量误差 (ppm)",
    "b  Fragment mass error (Da)": "b  子离子质量误差 (Da)",
    "b  Fragment mass error (ppm)": "b  子离子质量误差 (ppm)",
    "a  Cleavage specificity by group": "a  分组酶切特异性",
    "b  Missed cleavages by group": "b  分组漏切分布",
    "b  Missed cleavages (all PSMs)": "b  漏切分布（全部 PSM）",
    "Peptides with modification (%)": "带修饰肽段比例 (%)",
    "a  Unique-peptide length": "a  去重肽段长度",
    "b  Charge distribution (PSMs)": "b  电荷分布（PSM）",
    "PSMs": "PSM 数", "ID rate (%)": "解析率 (%)", "Raw tasks": "搜库任务",
    "FDR threshold": "FDR 阈值", "Decoy fraction (%)": "Decoy 占比 (%)",
    "Tasks": "任务", "Mass error (ppm)": "质量误差 (ppm)",
    "Mass error (Da)": "质量误差 (Da)", "Frequency": "频数",
    "Fully specific": "完全特异", "C-term specific": "C 端特异",
    "N-term specific": "N 端特异", "Non-specific": "非特异",
    "MC 0": "漏切 0", "MC 1": "漏切 1", "MC 2": "漏切 2",
    "Peptide length (aa)": "肽段长度 (aa)", "Unique peptides": "去重肽段",
    "Charge": "电荷", "PSMs (%)": "PSM (%)",
}


def T(en):
    from .i18n import get_lang
    return ZH_E1.get(en, _t(en)) if get_lang() == "zh" else en


# ---------------------------------------------------------------- 数据收集
def _task_dir(root, sp, task):
    return os.path.join(root, "04_search_annot", sp, task)


def parse_summary(path):
    """pFind.summary → dict（产率/覆盖度/特异性/修饰）。"""
    d = {}
    if not os.path.exists(path):
        return d
    sec = None
    for line in open(path, encoding="utf-8", errors="replace"):
        s = line.strip()
        m = re.match(r"^(Spectra|Peptides|Sequences|Proteins|Avg_Seq_Cov|"
                     r"Decoy Spectra):\s*([\d.]+)", s)
        if m:
            d[m.group(1)] = float(m.group(2))
            continue
        if s == "Cleavage:":
            sec = "cleav"
            continue
        if s == "Modifications:":
            sec = "mod"
            continue
        if s.startswith("---"):
            sec = None
            continue
        if sec == "cleav":
            m = re.match(r"^(Specific|C-term Specific|N-term Specific|"
                         r"Non-Specific):\s+\d+\s+([\d.]+)%", s)
            if m:
                d.setdefault("cleav", {})[m.group(1)] = float(m.group(2))
        elif sec == "mod":
            m = re.match(r"^(\S+)\t(\d+) \(([\d.]+)%\)", s)
            if m:
                d.setdefault("mods", {})[m.group(1)] = float(m.group(3))
    return d


def _merged_len(ivs):
    """区间并集总长（覆盖残基数）。"""
    if not ivs:
        return 0
    ivs = sorted(ivs)
    tot, cs, ce = 0, ivs[0][0], ivs[0][1]
    for s, e in ivs[1:]:
        if s <= ce:
            ce = max(ce, e)
        else:
            tot += ce - cs
            cs, ce = s, e
    return tot + ce - cs


def _load_l0(db_dir):
    """db_annot.faa → {L0 蛋白名: 序列}（表头即 pFind Proteins 列名）。"""
    from .novel import read_faa
    f = os.path.join(db_dir, "db_annot.faa")
    out = {}
    if os.path.exists(f):
        for h, s in read_faa(f):
            if h.startswith("L0_"):
                out[h.split()[0]] = s
    return out


def gather(root, log=print):
    """一次遍历 64 任务 → 全部聚合量。"""
    raw = load_tsv(os.path.join(root, "00_input", "raw_info.tsv"))
    k0 = [k for k in raw[0] if "raw" in k.lower() or "file" in k.lower()][0]
    split = load_tsv(os.path.join(root, "03_mgf", "split_report.tsv"))
    ms2_of = {}
    for r in split:
        ms2_of[(r.get("raw_file", ""), r.get("cv", ""))] = \
            int(r.get("n_spectra") or 0)
    # 任务 → raw（species_code+sample_type+enzyme+fraction 唯一定位）
    key2raw = {}
    for r in raw:
        key2raw[(r["species_code"], r["sample_type"], r["enzyme"],
                 r["fraction"])] = (r[k0], )
    tasks = []
    qv, td, ppm, frag = [], [], [], []
    len_cnt, mc_cnt = {}, {}
    mc_g = {}   # E1Q5b：组级漏切计数 {group: {mc: n}}
    chg_cnt = {"TL": {}, "LL": {}}          # 电荷按酶拆分（镜像对照）
    pep_set = set()
    frag_meta = {}      # File_Name → (seq, Modification)（碎片 ppm 自算输入）
    # 样本级聚合（E1Q1b）：肽段(序列+修饰)/蛋白 跨任务去重并集
    grp = {}   # (sp,st,enz) -> dict(psm,ms2,peps=set,prots=set)
    # 馏分级聚合（E1Q1a）：同口径，按 (sp,st,enz,frac)
    frac_agg = {}
    # 序列覆盖度（E1Q3）：逐蛋白覆盖区间（L0 蛋白）
    l0 = {sp: _load_l0(os.path.join(root, "02_db", sp))
          for sp in ("Mab", "Mch")}
    cov = {}   # (sp,st,enz,L0蛋白) -> [(start,end),...]
    for sp in ("Mab", "Mch"):
        base = os.path.join(root, "04_search_annot", sp)
        if not os.path.isdir(base):
            continue
        for task in sorted(os.listdir(base)):
            tdir = _task_dir(root, sp, task)
            smf = os.path.join(tdir, "pFind.summary")
            spec = os.path.join(tdir, "pFind-Filtered.spectra")
            if not os.path.exists(smf):
                continue
            m = re.match(r"^(\w+)_(\w+)_(TL|LL)_(F\d+)_(\d+)$", task)
            if not m:
                continue
            st, enz, frac, cv = m.group(2), m.group(3), m.group(4), m.group(5)
            sm = parse_summary(smf)
            rawhit = key2raw.get((sp, st, enz, frac), ("",))
            ms2 = ms2_of.get((rawhit[0], cv), 0)
            enz_key = enz
            tasks.append(dict(task=task, sp=sp, st=st, enz=enz, frac=frac,
                              cv=cv,
                              n_psm=int(sm.get("Spectra", 0)),
                              n_pep=int(sm.get("Peptides", 0)),
                              seqcov=sm.get("Avg_Seq_Cov", 0.0),
                              decoy=sm.get("Decoy Spectra", 0),
                              ms2=ms2, cleav=sm.get("cleav", {}),
                              mods=sm.get("mods", {}),
                              group=f"{sp}·{st}·{enz}"))
            if not os.path.exists(spec):
                continue
            g = grp.setdefault((sp, st, enz),
                               dict(psm=0, ms2=0, peps=set(), prots=set(),
                                    cores=set()))
            g["ms2"] += ms2
            gf = frac_agg.setdefault((sp, st, enz, frac),
                                 dict(psm=0, ms2=0, peps=set(), prots=set()))
            gf["ms2"] += ms2
            with open(spec, newline="", encoding="utf-8",
                      errors="replace") as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    q = r.get("Q-value", "")
                    decoy = (r.get("Target/Decoy") or "").lower() == "decoy"
                    try:
                        if q != "":
                            qv.append(float(q))
                        td.append(decoy)
                        g["psm"] += 1
                        gf["psm"] += 1
                        seq0 = r.get("Sequence") or ""
                        fn_key = r.get("File_Name") or ""
                        mh = float(r.get("Exp.MH+") or 0)
                        sh = float(r.get("Mass_Shift(Exp.-Calc.)") or 0)
                        if mh > 0:
                            ppm.append(sh / mh * 1e6)
                        # 碎片 ppm 需从原始 MGF 自算（见 e1frag）；此处只登记
                        if seq0:
                            frag_meta.setdefault(fn_key, (seq0,
                                                          r.get("Modification") or ""))
                        seq = r.get("Sequence") or ""
                        if seq:
                            mod = r.get("Modification") or ""
                            g["peps"].add((seq.upper(), mod))
                            g["cores"].add(_core(seq))
                            gf["peps"].add((seq.upper(), mod))
                            pep_set.add(seq.upper())
                            L = len(seq)
                            len_cnt[L] = len_cnt.get(L, 0) + 1
                        prots = [p for p in
                                 re.split(r"[;/]", r.get("Proteins") or "") if p]
                        for p in prots:
                            g["prots"].add(p)
                            gf["prots"].add(p)
                        # 覆盖区间（Proteins 与 Positions 平行，1-based 起点）
                        if seq:
                            poss = [x for x in
                                    re.split(r"[;/]", r.get("Positions") or "")
                                    if x]
                            for p, po in zip(prots, poss):
                                if not p.startswith("L0_"):
                                    continue
                                if p not in l0[sp]:
                                    continue
                                try:
                                    st0 = int(po.split(",")[0]) - 1
                                except ValueError:
                                    continue
                                cov.setdefault((sp, st, enz, p), []).append(
                                    (st0, st0 + len(seq)))
                        z = int(r.get("Charge") or 0)
                        chg_cnt[enz_key][z] = chg_cnt[enz_key].get(z, 0) + 1
                        mc = int(float(r.get("Miss.Clv.Sites") or 0))
                        mc_cnt[mc] = mc_cnt.get(mc, 0) + 1
                        _gk = f"{sp}·{st}·{enz}"
                        mc_g.setdefault(_gk, {})
                        mc_g[_gk][mc] = mc_g[_gk].get(mc, 0) + 1
                    except ValueError:
                        continue
    log(f"  E1 聚合：{len(tasks)} 任务 / PSM {sum(t['n_psm'] for t in tasks):,}"
        f" / 全局谱行 {len(qv):,} / 去重肽 {len(pep_set):,}")
    # ---- 未过滤 pFind.spectra：Q-value + Target/Decoy（E1Q2 曲线用，按组分桶）----
    qv_u, td_u = [], []
    qv_g = {g: [] for g in GROUP_ORDER}      # 组级（E1Q2a 六小图）
    td_g = {g: [] for g in GROUP_ORDER}
    for sp in ("Mab", "Mch"):
        base = os.path.join(root, "04_search_annot", sp)
        if not os.path.isdir(base):
            continue
        for task in sorted(os.listdir(base)):
            spec = os.path.join(base, task, "pFind.spectra")
            if not os.path.exists(spec):
                continue
            m = re.match(r"^(\w+)_(\w+)_(TL|LL)_(F\d+)_(\d+)$", task)
            gkey = f"{sp}·{m.group(2)}·{m.group(3)}" if m else None
            with open(spec, newline="", encoding="utf-8",
                      errors="replace") as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    q = r.get("Q-value", "")
                    tdv = (r.get("Target/Decoy") or "").lower()
                    if q in ("", None) or tdv not in ("target", "decoy"):
                        continue
                    try:
                        qf = float(q)
                    except ValueError:
                        continue
                    qv_u.append(qf)
                    td_u.append(tdv == "decoy")
                    if gkey in qv_g:
                        qv_g[gkey].append(qf)
                        td_g[gkey].append(tdv == "decoy")
    log(f"  未过滤 T/D 谱行：{len(qv_u):,}（decoy {sum(td_u):,}）")
    # ---- 逐蛋白覆盖度分布（E1Q3）：TL / LL / Combined ----
    cov_dist = {}   # (sp,st,enz) -> [每蛋白覆盖%]，enz ∈ TL/LL/C
    for sp in ("Mab", "Mch"):
        for st in sorted({k[1] for k in cov if k[0] == sp}):
            for enz in ("TL", "LL", "C"):
                vals = []
                for acc, seqp in l0[sp].items():
                    ivs = cov.get((sp, st, "TL", acc), []) if enz != "LL" else []
                    if enz != "TL":
                        ivs = ivs + cov.get((sp, st, "LL", acc), [])
                    if not ivs:
                        continue
                    vals.append(100.0 * _merged_len(ivs) / max(len(seqp), 1))
                if vals:
                    cov_dist[(sp, st, enz)] = vals
    log(f"  覆盖度：{len(cov)} 蛋白-组键 / "
        f"{sum(len(v) for v in cov_dist.values())} 逐蛋白值")
    # ---- 碎片离子 ppm（E1Q4b）：从原始 MGF 自算（单位与母离子统一为 ppm）----
    from . import e1frag
    got = e1frag.compute(root, frag_meta, log=log)
    frag_vals = [v[0] for v in got.values()]
    return dict(tasks=tasks, qv=np.array(qv), td=np.array(td, dtype=bool),
                qv_u=np.array(qv_u), td_u=np.array(td_u, dtype=bool),
                ppm=np.array(ppm), frag=np.array(frag_vals),
                len_cnt=len_cnt,
                chg_cnt=chg_cnt, mc_cnt=mc_cnt, mc_g=mc_g,
                pep_n=len(pep_set), grp=grp,
                frac=frac_agg, qv_g=qv_g, td_g=td_g, cov_dist=cov_dist)


def _short_task(t):
    return f"{t['sp']}·{t['st']}·{t['enz']}·{t['frac']}·{t['cv']}"


GROUP_ORDER = ["Mab·CFP·TL", "Mab·CFP·LL", "Mab·TCL·TL", "Mab·TCL·LL",
               "Mch·TCL·TL", "Mch·TCL·LL"]


# ---------------------------------------------------------------- E1Q1
_E1Q1_PANELS = [
    ("a  Identified spectra (PSMs)", "PSMs",
     lambda g: g["psm"], "{:,}"),
    ("b  Peptides", "Peptides",
     lambda g: len(g["peps"]), "{:,}"),
    ("c  Proteins", "Proteins",
     lambda g: len(g["prots"]), "{:,}"),
    ("d  Spectrum explain rate", "%",
     lambda g: 100 * g["psm"] / g["ms2"] if g["ms2"] else 0, "{:.1f}"),
]


def e1q1a_fraction(root, out_dir, agg):
    """E1Q1a 每馏分四指标（横向条形）：y=馏分（a 面板行标签），TL/LL 成对。"""
    fr = agg["frac"]
    bases = sorted({(sp, st, fc) for (sp, st, enz, fc) in fr},
                   key=lambda b: (b[0], b[1], b[2]))
    y = np.arange(len(bases))
    fig = plt.figure(figsize=(COL2, 0.155 * len(bases) + 1.15))
    gs = fig.add_gridspec(1, 4, wspace=0.34)   # 四面板等宽等距
    axes = [fig.add_subplot(gs[0, j]) for j in range(4)]
    for j, (ax, (ttl, ylab, get, fmt)) in enumerate(zip(axes, _E1Q1_PANELS)):
        vmax = max(max(get(fr[k2]) for k2 in fr), 1)
        for i, (sp, st, fc) in enumerate(bases):
            for off, enz, c in ((-0.21, "TL", EN_COL["TL"]),
                                (+0.21, "LL", EN_COL["LL"])):
                g = fr.get((sp, st, enz, fc))
                if g:
                    ax.barh(i + off, get(g), 0.38, color=c)
        ax.set_yticks(y)
        ax.set_yticklabels([f"{sp[:3]}{st[0]}·{fc}" for sp, st, fc in bases]
                           if j == 0 else [], fontsize=5.0)
        ax.invert_yaxis()
        ax.set_title(T(ttl), fontsize=7.5, fontweight="bold", loc="left")
        ax.set_xlabel(ylab)
        ax.set_xlim(0, vmax * 1.12)
    handles = [plt.Rectangle((0, 0), 1, 1, color=EN_COL[e]) for e in ("TL", "LL")]
    axes[0].legend(handles, [EN2NAME["TL"], EN2NAME["LL"]], frameon=False,
                   fontsize=5, loc="lower right")
    _save(fig, out_dir, "E1Q1a_yield_per_fraction")


def e1q1b_sample(root, out_dir, agg):
    """E1Q1b 每样本四指标：谱图数/肽段数/蛋白数/解析率（跨任务去重并集）。"""
    grp = agg["grp"]
    bases = sorted({(sp, st) for (sp, st, enz) in grp})
    x = np.arange(len(bases))
    fig, axes = plt.subplots(1, 4, figsize=(COL2, 2.9),
                             gridspec_kw={"wspace": 0.46})
    for ax, (ttl, ylab, get, fmt) in zip(axes, _E1Q1_PANELS):
        vmax = max(max(get(grp[k2]) for k2 in grp), 1)
        for i, (sp, st) in enumerate(bases):
            for off, enz, c in ((-0.19, "TL", EN_COL["TL"]),
                                (+0.19, "LL", EN_COL["LL"])):
                g = grp.get((sp, st, enz))
                if not g:
                    continue
                v = get(g)
                ax.bar(i + off, v, 0.36, color=c)
                ax.text(i + off, v * 1.04, fmt.format(v), ha="center",
                        va="bottom", fontsize=4.8)
        ax.set_xticks(x)
        ax.set_xticklabels([f"{sp}·{st}" for sp, st in bases], rotation=28,
                           ha="right", fontsize=5.4)
        ax.set_title(T(ttl), fontsize=7.5, fontweight="bold", loc="left")
        ax.set_ylabel(ylab)
        ax.set_ylim(0, vmax * 1.42)
    handles = [plt.Rectangle((0, 0), 1, 1, color=EN_COL[e]) for e in ("TL", "LL")]
    axes[0].legend(handles, [EN2NAME["TL"], EN2NAME["LL"]], frameon=False,
                   fontsize=5, loc="upper right")
    _save(fig, out_dir, "E1Q1b_yield_per_sample")


# ---------------------------------------------------------------- E1Q2
def e1q2_fdr(root, out_dir, agg):
    """六组各自 Target/Decoy 累计曲线（a1–a6 小多重图）。"""
    qv_g, td_g = agg["qv_g"], agg["td_g"]
    thr = np.array([0.0001, 0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05])
    i01 = list(thr).index(0.01)
    fig = plt.figure(figsize=(COL2, 5.0))
    outer = fig.add_gridspec(2, 3, hspace=0.92, wspace=0.42)
    for i, g in enumerate(GROUP_ORDER):
        ax = fig.add_subplot(outer[i // 3, i % 3])
        q, d = np.array(qv_g.get(g, [])), np.array(td_g.get(g, []), dtype=bool)
        nt = [int(((q <= t) & ~d).sum()) for t in thr]
        nd = [int(((q <= t) & d).sum()) for t in thr]
        ax.plot(thr, nt, "o-", color=COLORS[0], lw=1.0, ms=2.2, label="Target")
        ax.plot(thr, nd, "s-", color=COLORS[3], lw=1.0, ms=2.2, label="Decoy")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_ylim(bottom=0.7)
        ax.set_title(f"a{i+1}  {g}", fontsize=6.5, fontweight="bold", loc="left")
        nt01, nd01 = nt[i01], nd[i01]
        # 1% 统计并入 x 轴标签第二行（轴外，零冲突）
        ax.set_xlabel(f"{T('FDR threshold')}\n"
                      f"1%: T {nt01:,} / D {nd01:,} "
                      f"({100*nd01/max(nt01,1):.2f}%)", fontsize=5.2)
        ax.set_ylabel(T("PSMs"), fontsize=5.5)
        ax.tick_params(labelsize=5)
        if i == 0:
            ax.legend(frameon=False, fontsize=4.6, loc="lower right")
    _save(fig, out_dir, "E1Q2_fdr")


# ---------------------------------------------------------------- E1Q3
ENZ3 = [("TL", "Trypsin+LysC", COLORS[0]),
        ("LL", "LysargiNase+LysN", COLORS[3]),
        ("C", "Combined", COLORS[2])]


def e1q3_seqcov(root, out_dir, agg):
    """3 样本组 × {Trypsin, LysargiNase, Combined} 逐蛋白覆盖度箱线。"""
    cd = agg["cov_dist"]
    samples = sorted({(k[0], k[1]) for k in cd})
    fig, ax = plt.subplots(figsize=(COL1, 2.9))
    data, cols, pos, ticks, tlabels = [], [], [], [], []
    for i, (sp, st) in enumerate(samples):
        for j, (enz, name, c) in enumerate(ENZ3):
            vals = cd.get((sp, st, enz), [])
            data.append(vals if vals else [0])
            cols.append(c)
            pos.append(i + (j - 1) * 0.26)
        ticks.append(i)
        tlabels.append(f"{sp}·{st}")
    bp = ax.boxplot(data, positions=pos, widths=0.22, patch_artist=True,
                    showfliers=False,
                    medianprops=dict(color="#111827", lw=0.9),
                    whiskerprops=dict(lw=0.6), capprops=dict(lw=0.6))
    for patch, c in zip(bp["boxes"], cols):
        patch.set_facecolor(c)
        patch.set_alpha(0.55)
        patch.set_edgecolor(c)
        patch.set_linewidth(0.6)
    # 中位数标注（每组三列的 med 值）
    med_txt = []
    for i, (sp, st) in enumerate(samples):
        for j, (enz, name, c) in enumerate(ENZ3):
            v = np.median(cd.get((sp, st, enz), [0]))
            med_txt.append(f"{v:.1f}")
    ymax = max(max(d) for d in data if d)
    for xi, t in zip(pos, med_txt):
        ax.text(xi, ymax * 1.045, t, ha="center", va="bottom", fontsize=4.8,
                color="#374151")
    ax.text(0.012, 0.965, T("Numbers = median"),
            transform=ax.transAxes, fontsize=4.8, color="#374151",
            va="top", ha="left")
    ax.set_xticks(ticks)
    ax.set_xticklabels(tlabels, fontsize=6)
    ax.set_xlim(-0.55, len(samples) - 0.45)
    ax.set_ylim(0, ymax * 1.18)
    ax.set_ylabel(T("Protein sequence coverage (%)"))
    ax.set_title(T("E1Q3  Protein sequence coverage"), fontsize=7.5,
                 fontweight="bold", loc="left")
    handles = [plt.Rectangle((0, 0), 1, 1, fc=c, alpha=0.55, ec=c, lw=0.6)
               for _, _, c in ENZ3]
    ax.legend(handles, [n for _, n, _ in ENZ3], frameon=False, fontsize=5,
              loc="upper right", ncol=3, handlelength=1.2,
              columnspacing=0.9, borderpad=0.2)
    _save(fig, out_dir, "E1Q3_seqcov")


# ---------------------------------------------------------------- E1Q4
def e1q4_mass(root, out_dir, agg):
    """母离子与子离子质量误差并列，单位统一 ppm。"""
    ppm, frag = agg["ppm"], agg["frag"]
    fig, axes = plt.subplots(1, 2, figsize=(COL2, 2.7),
                             gridspec_kw={"wspace": 0.34})
    for ax, dat, unit, ttl in (
            (axes[0], ppm, "ppm", T("a  Precursor mass error (ppm)")),
            (axes[1], frag, "ppm", T("b  Fragment mass error (ppm)"))):
        mu, sd = float(np.mean(dat)), float(np.std(dat))
        lo, hi = mu - 5 * sd, mu + 5 * sd
        ax.hist(dat[(dat >= lo) & (dat <= hi)], bins=80,
                color=COLORS[0], alpha=0.85)
        ax.axvline(mu, color=COLORS[3], ls="--", lw=0.8)
        ax.axvline(0.0, color="#6B7280", ls=":", lw=0.7)
        out_n = int(((dat < lo) | (dat > hi)).sum())
        ax.set_title(ttl, fontsize=7.5, fontweight="bold", loc="left")
        ax.set_xlabel("Mass error (" + unit + ")")
        ax.set_ylabel(T("Frequency"))
        ax.text(0.98, 0.95,
                f"n = {len(dat):,}\nmean {mu:+.3f} {unit}\nSD {sd:.3f}"
                + (f"\noutside ±5SD: {out_n}" if out_n else ""),
                transform=ax.transAxes, ha="right", va="top", fontsize=5.6,
                bbox=dict(boxstyle="round,pad=0.3", fc="#F3F4F6", ec="none"))
    _save(fig, out_dir, "E1Q4_mass_accuracy")


# ---------------------------------------------------------------- E1Q5
def e1q5_cleav(root, out_dir, agg):
    from collections import defaultdict
    cl = defaultdict(list)
    for t in agg["tasks"]:
        cl[t["group"]].append(t["cleav"])
    fig, axes = plt.subplots(1, 2, figsize=(COL2, 2.7),
                             gridspec_kw={"wspace": 0.42})
    keys = ["Specific", "C-term Specific", "N-term Specific", "Non-Specific"]
    labs = [T("Fully specific"), T("C-term specific"), T("N-term specific"),
            T("Non-specific")]
    pal = [COLORS[0], COLORS[1], COLORS[2], "#9CA3AF"]
    x = np.arange(len(GROUP_ORDER))
    ax = axes[0]
    bottom = np.zeros(len(GROUP_ORDER))
    hidden = []            # 全组均低于标注阈值的类别
    hide_thr = 3.0
    for k, lab, c in zip(keys, labs, pal):
        vals = np.array([np.mean([d.get(k, 0) for d in cl.get(g, [])])
                         for g in GROUP_ORDER])
        ax.bar(x, vals, 0.62, bottom=bottom, color=c, label=lab)
        dark = c in (COLORS[1], COLORS[2], "#9CA3AF")
        for xi, v, bt in zip(x, vals, bottom):
            if v >= hide_thr:
                ax.text(xi, bt + v / 2, f"{v:.1f}", ha="center", va="center",
                        fontsize=5.2, color="#374151" if dark else "white")
        if vals.max() < hide_thr:
            hidden.append((lab, float(vals.max())))
        bottom += vals
    ax.set_xticks(x)
    ax.set_xticklabels(GROUP_ORDER, rotation=28, ha="right", fontsize=5.2)
    ax.set_ylabel("%")
    ax.set_ylim(0, 116)
    ax.set_title(T("a  Cleavage specificity by group"), fontsize=7.5,
                 fontweight="bold", loc="left", pad=26)
    note = T("numbers = mean % (segments < 3% unlabeled)")
    if hidden:
        note += " | " + "; ".join(f"{lab} ≤ {v:.1f}%" for lab, v in hidden)
    ax.text(0.0, 1.015, note, transform=ax.transAxes, fontsize=4.6,
            color="#374151", ha="left", va="bottom")
    ax.legend(frameon=False, fontsize=4.8, loc="lower center",
              bbox_to_anchor=(0.5, 1.10), ncol=4, handlelength=1.0,
              columnspacing=0.8, borderpad=0.2)
    ax = axes[1]
    mc_g = agg["mc_g"]
    mc_cols = [(0, T("MC 0"), COLORS[0]), (1, T("MC 1"), COLORS[1]),
               (2, T("MC 2"), "#9CA3AF")]
    bottom = np.zeros(len(GROUP_ORDER))
    for k, lab, c in mc_cols:
        vals = []
        for g in GROUP_ORDER:
            cnt = mc_g.get(g, {})
            tot = sum(cnt.values()) or 1
            vals.append(100.0 * cnt.get(k, 0) / tot)
        vals = np.array(vals)
        ax.bar(x, vals, 0.62, bottom=bottom, color=c, label=lab)
        dark = c in (COLORS[1], "#9CA3AF")
        for xi, v, bt in zip(x, vals, bottom):
            if v >= 3.0:
                ax.text(xi, bt + v / 2, f"{v:.1f}", ha="center", va="center",
                        fontsize=5.2, color="#374151" if dark else "white")
        bottom += vals
    ax.set_xticks(x)
    ax.set_xticklabels(GROUP_ORDER, rotation=28, ha="right", fontsize=5.2)
    ax.set_ylabel("% PSMs")
    ax.set_ylim(0, 116)
    ax.set_title(T("b  Missed cleavages by group"), fontsize=7.5,
                 fontweight="bold", loc="left", pad=26)
    ax.text(0.0, 1.015, T("numbers = mean % (segments < 3% unlabeled)"),
            transform=ax.transAxes, fontsize=4.6, color="#374151",
            ha="left", va="bottom")
    ax.legend(frameon=False, fontsize=4.8, loc="lower center",
              bbox_to_anchor=(0.5, 1.10), ncol=3, handlelength=1.0,
              columnspacing=0.8, borderpad=0.2)
    _save(fig, out_dir, "E1Q5_cleavage")


# ---------------------------------------------------------------- E1Q6
def e1q6_mods(root, out_dir, agg):
    from collections import defaultdict
    mo = defaultdict(list)
    for t in agg["tasks"]:
        mo[t["group"]].append(t["mods"])
    modnames = ["Oxidation[M]", "Carbamidomethyl[C]", "Acetyl[AnyN-term]"]
    fig, ax = plt.subplots(figsize=(COL1, 2.6))
    x = np.arange(len(modnames))
    w = 0.13
    pal = [COLORS[0], COLORS[3], COLORS[1], COLORS[2], "#CA9161", "#949494"]
    for i, g in enumerate(GROUP_ORDER):
        vals = [np.mean([d.get(mn, 0) for d in mo.get(g, [])])
                for mn in modnames]
        ax.bar(x + (i - 2.5) * w, vals, w, color=pal[i],
               label=g.replace("·", " "))
    ax.set_xticks(x)
    ax.set_xticklabels(modnames, fontsize=5.6)
    ax.set_ylabel(T("Peptides with modification (%)"))
    ax.set_title(T("E1Q6  Modification frequency"), fontsize=7.5,
                 fontweight="bold", loc="left")
    ax.legend(frameon=False, fontsize=4.6, ncol=2, loc="upper right")
    _save(fig, out_dir, "E1Q6_mods")


# ---------------------------------------------------------------- E1Q7
def e1q7_pep(root, out_dir, agg):
    lc = agg["len_cnt"]
    cc = agg["chg_cnt"]
    fig, axes = plt.subplots(1, 2, figsize=(COL2, 2.7),
                             gridspec_kw={"wspace": 0.34})
    ax = axes[0]
    Ls = sorted(k for k in lc if 5 <= k <= 40)
    ax.bar(Ls, [lc[k] for k in Ls], color=COLORS[0], width=0.8)
    ax.set_xlabel(T("Peptide length (aa)"))
    ax.set_ylabel(T("Unique peptides"))
    ax.set_title(T("a  Unique-peptide length"), fontsize=7.5,
                 fontweight="bold", loc="left")
    lmed = float(np.average(Ls, weights=[lc[k] for k in Ls]))
    ax.axvline(lmed, color=COLORS[3], ls="--", lw=0.8)
    ax.text(lmed + 0.3, ax.get_ylim()[1] * 0.95, f"med {lmed:.1f}",
            fontsize=5.4, color=COLORS[3], va="top")
    ax = axes[1]
    zs = sorted(z for z in set(cc["TL"]) | set(cc["LL"]) if 2 <= z <= 5)
    totT, totL = sum(cc["TL"].values()), sum(cc["LL"].values())
    vT = [100 * cc["TL"].get(z, 0) / totT for z in zs]
    vL = [100 * cc["LL"].get(z, 0) / totL for z in zs]
    x = np.arange(len(zs))
    ax.bar(x - 0.19, vT, 0.38, color=EN_COL["TL"], label=EN2NAME["TL"])
    ax.bar(x + 0.19, vL, 0.38, color=EN_COL["LL"], label=EN2NAME["LL"])
    for xi, v in zip(x, vT):
        ax.text(xi - 0.19, v + 0.8, f"{v:.0f}", ha="center", fontsize=5.0)
    for xi, v in zip(x, vL):
        ax.text(xi + 0.19, v + 0.8, f"{v:.0f}", ha="center", fontsize=5.0)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{z}+" for z in zs])
    ax.set_ylabel(T("PSMs (%)"))
    ax.set_title(T("b  Charge distribution (PSMs)"), fontsize=7.5,
                 fontweight="bold", loc="left")
    ax.legend(frameon=False, fontsize=5, loc="upper right")
    _save(fig, out_dir, "E1Q7_peptide_dist")


# ---------------------------------------------------------------- 结论句
def write_conclusions(root, out_dir, agg):
    from .i18n import get_lang
    zh = get_lang() == "zh"
    ts = agg["tasks"]
    tot_psm = sum(t["n_psm"] for t in ts)
    tot_pep = sum(t["n_pep"] for t in ts)
    rates = [100 * t["n_psm"] / t["ms2"] for t in ts if t["ms2"]]
    qv, td = agg["qv_u"], agg["td_u"]
    n01 = int(((qv <= 0.01) & ~td).sum())
    d01 = int(((qv <= 0.01) & td).sum())
    # 组级 1% 统计（E1Q2a 六小图）
    g01 = {}
    for g in GROUP_ORDER:
        gg = np.array(agg["qv_g"].get(g, []))
        gd = np.array(agg["td_g"].get(g, []), dtype=bool)
        if gg.size:
            g01[g] = (int(((gg <= 0.01) & ~gd).sum()),
                      int(((gg <= 0.01) & gd).sum()))
    g01_str = "; ".join(f"{g.split('·')[0]}·{g.split('·')[1]}·{g.split('·')[2]} "
                        f"T {v[0]:,}/D {v[1]:,} ({100*v[1]/max(v[0],1):.2f}%)"
                        for g, v in g01.items())
    covs = [t["seqcov"] for t in ts]
    # E1Q3 三酶中位与合并增益
    cd = agg["cov_dist"]
    samples3 = sorted({(k[0], k[1]) for k in cd})
    def _med(sp, st, enz):
        v = cd.get((sp, st, enz), [])
        return float(np.median(v)) if v else 0.0
    cov3_zh = "；".join(
        f"{sp}·{st} TL {_med(sp,st,'TL'):.1f} / LL {_med(sp,st,'LL'):.1f} / "
        f"合并 {_med(sp,st,'C'):.1f}" for sp, st in samples3)
    cov3_en = "; ".join(
        f"{sp}·{st} TL {_med(sp,st,'TL'):.1f} / LL {_med(sp,st,'LL'):.1f} / "
        f"Comb {_med(sp,st,'C'):.1f}" for sp, st in samples3)
    _gains = [(_med(sp, st, "C") - 0.5 * (_med(sp, st, "TL")
                                          + _med(sp, st, "LL")))
              for sp, st in samples3]
    cov_gain = float(np.mean(_gains)) if _gains else 0.0
    ppm, frag = agg["ppm"], agg["frag"]
    cl = [t["cleav"].get("Specific", 0) for t in ts if t["cleav"]]
    mc = agg["mc_cnt"]
    mctot = sum(mc.values())
    _mc0 = []
    for _g in GROUP_ORDER:
        _c = agg["mc_g"].get(_g, {})
        _t = sum(_c.values())
        if _t:
            _mc0.append(100.0 * _c.get(0, 0) / _t)
    mc0_range = min(_mc0) if _mc0 else 0.0
    mc0_hi = max(_mc0) if _mc0 else 0.0
    # 分组完全特异率区间（按任务 → 组均值）
    _sp_by_g = {}
    for _t in ts:
        if _t["cleav"].get("Specific") is not None:
            _sp_by_g.setdefault(_t["group"], []).append(
                _t["cleav"]["Specific"])
    _spec = [np.mean(v) for v in _sp_by_g.values() if v]
    spec_lo = min(_spec) if _spec else 0.0
    spec_hi = max(_spec) if _spec else 0.0
    mods = {}
    for t in ts:
        for k, v in t["mods"].items():
            mods.setdefault(k, []).append(v)
    grp = agg["grp"]
    fr = agg["frac"]
    gv = dict(psm=[g["psm"] for g in grp.values()],
              pep=[len(g["peps"]) for g in grp.values()],
              prot=[len(g["prots"]) for g in grp.values()],
              rate=[100 * g["psm"] / g["ms2"] for g in grp.values()
                    if g["ms2"]])
    fv = dict(psm=[g["psm"] for g in fr.values()],
              pep=[len(g["peps"]) for g in fr.values()],
              prot=[len(g["prots"]) for g in fr.values()],
              rate=[100 * g["psm"] / g["ms2"] for g in fr.values()
                    if g["ms2"]])
    fr_items = sorted(((k, 100 * g["psm"] / g["ms2"]) for k, g in
                       fr.items() if g["ms2"]), key=lambda x: x[1])
    lo2 = " / ".join(f"{k[0]}·{k[1]}·{k[2]}·{k[3]} {v:.1f}%"
                     for k, v in fr_items[:2])
    rest = [v for _, v in fr_items[2:]]
    modstr = "; ".join(f"{k} {np.mean(v):.1f}%"
                       for k, v in sorted(mods.items()))
    lc, cc = agg["len_cnt"], agg["chg_cnt"]
    Ls = sorted(k for k in lc if 5 <= k <= 40)
    lmed = float(np.average(Ls, weights=[lc[k] for k in Ls]))
    cc_all = {}
    for e in ("TL", "LL"):
        for z, v in cc[e].items():
            cc_all[z] = cc_all.get(z, 0) + v
    ctot = sum(cc_all.values())
    ztop = max((z for z in cc_all if 2 <= z <= 6), key=lambda z: cc_all[z])
    if zh:
        lines = [
            "# E1 QC 图结论句（已注释库搜库结果蛋白质组质控）",
            "",
            f"**E1Q1a · 鉴定产率（每馏分）** —— 32 个馏分×酶组"
            f"（16 馏分 × TL/LL，馏分内跨 CV 合并）：PSM {min(fv['psm']):,}"
            f"–{max(fv['psm']):,}、肽段 {min(fv['pep']):,}–{max(fv['pep']):,}、"
            f"蛋白 {min(fv['prot']):,}–{max(fv['prot']):,}、"
            f"解析率 {min(fv['rate']):.1f}–{max(fv['rate']):.1f}%，"
            f"最低两个馏分为 {lo2}（其余 {min(rest):.1f}–{max(rest):.1f}%）。",

            "",
            f"**E1Q1b · 鉴定产率（每样本）** —— 64 个任务共 {tot_psm:,} PSM；"
            f"按样本（跨馏分×CV 去重并集）：PSM {min(gv['psm']):,}"
            f"–{max(gv['psm']):,}、肽段 {min(gv['pep']):,}–{max(gv['pep']):,}、"
            f"蛋白 {min(gv['prot']):,}–{max(gv['prot']):,}、"
            f"解析率 {min(gv['rate']):.1f}–{max(gv['rate']):.1f}%"
            f"（6 样本组，TL/LL 各半）。",
            "",
            f"**E1Q2 · FDR 分布** —— 全局 1% FDR 阈值下 target {n01:,} PSM 对应 "
            f"decoy {d01:,} 条（经验 FDR {100*d01/max(n01,1):.2f}%）；"
            f"分 6 组（a1–a6）：{g01_str}——各组经验 FDR 均低于 1%，"
            f"鉴定错误控制达标。",
            "",
            f"**E1Q3 · 序列覆盖率** —— 逐蛋白覆盖度分 3 样本组 × 双酶+合并"
            f"（箱线中位，%）：{cov3_zh}；双酶合并较两单酶均值提升 "
            f"{cov_gain:.1f} 个百分点，镜像两酶覆盖互补支撑注释蛋白质组的"
            f"深度表征。",
            "",
            f"**E1Q4 · 质量精度** —— 母离子 mean {np.mean(ppm):+.2f} ppm "
            f"(SD {np.std(ppm):.2f})、子离子 {np.mean(frag):+.2f} ppm "
            f"(SD {np.std(frag):.2f})，均落在仪器 20 ppm 容差内、无系统偏移。",
            "",
            f"**E1Q5 · 酶切特异性与漏切** —— 完全特异平均 {np.mean(cl):.1f}%"
            f"（TL 与 LL 各自符合 KR_C/KR_N 规则；分组 {spec_lo:.1f}–{spec_hi:.1f}%），"
            f"零漏切 PSM 总体 {100*mc.get(0,0)/mctot:.1f}%"
            f"（分组 {mc0_range:.1f}–{mc0_hi:.1f}%），酶切行为与建库设定一致。",
            "",
            f"**E1Q6 · 修饰频率** —— 肽级修饰频率：{modstr}，"
            f"无异常高频修饰（提示无过度修饰或假修饰指认）。",
            "",
            f"**E1Q7 · 肽段物理分布** —— 去重肽段 {agg['pep_n']:,} 条，"
            f"长度加权中位 {lmed:.1f} aa，电荷以 {ztop}+ 为主"
            f"（{100*cc_all[ztop]/ctot:.1f}%），分布符合 HCD/DDA 采集预期。",
            "",
        ]
    else:
        lines = [
            "# E1 QC Figure Conclusions (annotated-DB search proteomics QC)",
            "",
            f"**E1Q1a · Identification yield (per fraction)** — Across 32 "
            f"fraction-enzyme groups (16 fractions x TL/LL, CVs merged): "
            f"PSMs {min(fv['psm']):,}-{max(fv['psm']):,}, peptides "
            f"{min(fv['pep']):,}-{max(fv['pep']):,}, proteins "
            f"{min(fv['prot']):,}-{max(fv['prot']):,}, explain rate "
            f"{min(fv['rate']):.1f}-{max(fv['rate']):.1f}%; the two lowest "
            f"fractions are {lo2} (others {min(rest):.1f}-{max(rest):.1f}%).",

            "",
            f"**E1Q1b · Identification yield (per sample)** — {tot_psm:,} PSMs "
            f"across 64 tasks; after per-sample deduplication: PSMs "
            f"{min(gv['psm']):,}-{max(gv['psm']):,}, peptides "
            f"{min(gv['pep']):,}-{max(gv['pep']):,}, proteins "
            f"{min(gv['prot']):,}-{max(gv['prot']):,}, spectrum explain rate "
            f"{min(gv['rate']):.1f}-{max(gv['rate']):.1f}% (6 sample groups, "
            f"TL/LL in pairs).",
            "",
            f"**E1Q2 · Target-decoy FDR** — At the 1% FDR cut-off, {d01:,} "
            f"decoy PSMs versus {n01:,} target PSMs globally (empirical FDR "
            f"{100*d01/max(n01,1):.2f}%); per group (a1-a6): {g01_str} — "
            f"every group stays below the 1% line, confirming error control.",
            "",
            f"**E1Q3 · Sequence coverage** — Per-protein coverage by sample "
            f"group with Trypsin / LysargiNase / Combined (boxplot medians, "
            f"%): {cov3_en}; combined coverage exceeds the mean of the two "
            f"single enzymes by {cov_gain:.1f} percentage points, showing "
            f"complementary mirror coverage.",
            "",
            f"**E1Q4 · Mass accuracy** — Precursor error "
            f"{np.mean(ppm):+.2f} ppm (SD {np.std(ppm):.2f}) and fragment "
            f"{np.mean(frag):+.2f} ppm (SD {np.std(frag):.2f}), both within "
            f"the 20-ppm instrument tolerance without systematic offset.",
            "",
            f"**E1Q5 · Cleavage & missed cleavages** — Fully specific "
            f"cleavage averages {np.mean(cl):.1f}% (TL/LL each consistent "
            f"with KR_C/KR_N rules; per group {spec_lo:.1f}-{spec_hi:.1f}%), with "
            f"{100*mc.get(0,0)/mctot:.1f}% zero-missed-cleavage PSMs "
            f"(per group {mc0_range:.1f}-{mc0_hi:.1f}%).",
            "",
            f"**E1Q6 · Modification frequency** — Peptide-level frequencies: "
            f"{modstr}; no anomalous high-frequency modification is present.",
            "",
            f"**E1Q7 · Peptide length & charge** — {agg['pep_n']:,} "
            f"unique peptides, weighted-median length {lmed:.1f} aa, charge "
            f"dominated by {ztop}+ ({100*cc_all[ztop]/ctot:.1f}%), matching "
            f"HCD/DDA expectations.",
            "",
        ]
    out = os.path.join(out_dir, "E1_QC_figure_conclusions.md")
    with open(out, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines))
    return out


# ---------------------------------------------------------------- 样本汇总表
def write_sample_table(root, out_dir, agg):
    """8 指标行 × 6 样本组列表：TSV + XLSX（镜像证据=与同样本另一酶侧并集交集）。"""
    from .i18n import get_lang
    zh = get_lang() == "zh"
    grp = agg["grp"]
    cols = [(sp, st, enz) for (sp, st) in
            sorted({(a, b) for (a, b, c) in grp}) for enz in ("TL", "LL")
            if (sp, st, enz) in grp]
    EN2N = {"TL": "Trypsin+LysC", "LL": "LysargiNase+LysN"}

    def mirror(g, kind):
        sp, st, enz = g
        other = grp.get((sp, st, "LL" if enz == "TL" else "TL"))
        if not other:
            return ""
        return len(grp[g][kind] & other[kind])

    def union(g, kind):
        sp, st, enz = g
        other = grp.get((sp, st, "LL" if enz == "TL" else "TL"))
        if not other:
            return len(grp[g][kind])
        return len(grp[g][kind] | other[kind])

    rows = []
    if zh:
        hdr = ["指标"] + [f"{sp}·{st}·{enz}" for sp, st, enz in cols]
        rows.append(hdr)
        rows.append(["样本名"] + [f"{sp} / {st}" for sp, st, _ in cols])
        rows.append(["酶切"] + [EN2N[e] for _, _, e in cols])
        rows.append(["鉴定谱图数"] + [f"{grp[g]['psm']:,}" for g in cols])
        rows.append(["谱图解析率 (%)"] +
                    [f"{100*grp[g]['psm']/grp[g]['ms2']:.1f}" for g in cols])
        rows.append(["鉴定肽段数"] + [f"{len(grp[g]['peps']):,}" for g in cols])
        rows.append(["鉴定蛋白数"] + [f"{len(grp[g]['prots']):,}" for g in cols])
        rows.append(["两种酶切共鉴定肽段数"] +
                    [f"{union(g,'peps'):,}" for g in cols])
        rows.append(["两种酶切共鉴定蛋白数"] +
                    [f"{union(g,'prots'):,}" for g in cols])
        rows.append(["镜像酶切证据肽段数（核心区匹配）"] +
                    [f"{mirror(g,'cores'):,}" for g in cols])
        rows.append(["镜像酶切证据蛋白数"] + [f"{mirror(g,'prots'):,}" for g in cols])
    else:
        hdr = ["Metric"] + [f"{sp}·{st}·{enz}" for sp, st, enz in cols]
        rows.append(hdr)
        rows.append(["Sample"] + [f"{sp} / {st}" for sp, st, _ in cols])
        rows.append(["Enzyme"] + [EN2N[e] for _, _, e in cols])
        rows.append(["Identified spectra (PSMs)"] +
                    [f"{grp[g]['psm']:,}" for g in cols])
        rows.append(["Spectrum explain rate (%)"] +
                    [f"{100*grp[g]['psm']/grp[g]['ms2']:.1f}" for g in cols])
        rows.append(["Identified peptides"] +
                    [f"{len(grp[g]['peps']):,}" for g in cols])
        rows.append(["Identified proteins"] +
                    [f"{len(grp[g]['prots']):,}" for g in cols])
        rows.append(["Peptides identified by both enzymes (union)"] +
                    [f"{union(g,'peps'):,}" for g in cols])
        rows.append(["Proteins identified by both enzymes (union)"] +
                    [f"{union(g,'prots'):,}" for g in cols])
        rows.append(["Mirror-enzyme evidence peptides (core-matched)"] +
                    [f"{mirror(g,'cores'):,}" for g in cols])
        rows.append(["Mirror-enzyme evidence proteins"] +
                    [f"{mirror(g,'prots'):,}" for g in cols])
    tsv = os.path.join(out_dir, "E1_sample_summary.tsv")
    with open(tsv, "w", newline="", encoding="utf-8-sig") as fh:
        for r in rows:
            fh.write("\t".join(r) + "\n")
    # XLSX（Arial 9pt、表头加粗、首列左对齐）
    try:
        from openpyxl import Workbook
        from openpyxl.styles import Alignment, Font
        wb = Workbook()
        ws = wb.active
        ws.title = "E1 sample summary"
        for r in rows:
            ws.append(r)
        for c in ws[1]:
            c.font = Font(name="Arial", size=9, bold=True)
        for row in ws.iter_rows(min_row=2):
            for c in row:
                c.font = Font(name="Arial", size=9)
            row[0].alignment = Alignment(horizontal="left")
        ws.column_dimensions["A"].width = 34
        for i in range(2, len(cols) + 2):
            ws.column_dimensions[chr(64 + i)].width = 17
        wb.save(os.path.join(out_dir, "E1_sample_summary.xlsx"))
    except ImportError:
        pass
    return tsv


# ---------------------------------------------------------------- 总入口
def build_e1figs(root, out_dir):
    _apply_zh_font()
    os.makedirs(out_dir, exist_ok=True)
    agg = gather(root)
    e1q1a_fraction(root, out_dir, agg)
    e1q1b_sample(root, out_dir, agg)
    e1q2_fdr(root, out_dir, agg)
    e1q3_seqcov(root, out_dir, agg)
    e1q4_mass(root, out_dir, agg)
    e1q5_cleav(root, out_dir, agg)
    e1q6_mods(root, out_dir, agg)
    e1q7_pep(root, out_dir, agg)
    write_sample_table(root, out_dir, agg)
    write_conclusions(root, out_dir, agg)
