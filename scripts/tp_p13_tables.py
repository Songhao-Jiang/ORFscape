# -*- coding: utf-8 -*-
"""TP 结果四件套（P13）：分类统计 / 质控分布 / 事件统计 / PDF 输入。

口径（用户 2026-11-11 定）：
  保留新肽段 = P7 `final_keep=KEEP` 的唯一肽（230 条 / 369 行）
  排除新肽段 = 其余新肽（690 − 230 = 460 条）
  已注释肽段 = E1（db_annot 搜库）里命中注释 CDS（L0_*）的肽
酶切：TL = Trypsin，LL = LysargiNase。

用法
----
    python scripts/tp_p13_tables.py --workdir <orfscape 项目> --out <输出目录>
"""
import argparse
import csv
import math
import os
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pandas as pd

from orfscape.core.novel import iln, load_spectra, read_faa

# ---------------------------------------------------------------- 常数
AA_MASS = {
    "G": 57.02146, "A": 71.03711, "S": 87.03203, "P": 97.05276, "V": 99.06841,
    "T": 101.04768, "C": 103.00919, "L": 113.08406, "I": 113.08406,
    "N": 114.04293, "D": 115.02694, "Q": 128.05858, "K": 128.09496,
    "E": 129.04259, "M": 131.04049, "H": 137.05891, "F": 147.06841,
    "R": 156.10111, "Y": 163.06333, "W": 186.07931,
}
PROTON = 1.0072765
H2O = 18.0105646
TOL = 20.0            # ppm
ENZ_NAME = {"TL": "Trypsin", "LL": "LysargiNase"}
ENZ_EN = {"TL": "Trypsin", "LL": "LysargiNase", "ALL": "Combined"}
CAT_EN = {"已注释": "Annotated", "排除新肽段": "Filtered", "保留新肽段": "Obtained"}


def b_y_ions(seq, z):
    """理论 b/y 离子 [(mz, label)]，电荷 1..z-1。"""
    out = []
    n = len(seq)
    pre, suf = 0.0, 0.0
    mass = [AA_MASS.get(a, 0.0) for a in seq]
    total = sum(mass)
    for i in range(n - 1):
        pre += mass[i]
        suf = total - pre
        for c in range(1, max(2, z)):
            b = (pre + PROTON * c) / c
            y = (suf + H2O + PROTON * c) / c
            out.append((b, "b%d^%d" % (i + 1, c)))
            out.append((y, "y%d^%d" % (n - i - 1, c)))
    return out


def read_mgf_peaks(paths, wanted):
    """流式读 MGF，只收 wanted 里的 (stem, scan)。"""
    out = {}
    for p in paths:
        title = None
        buf = []
        with open(p, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                line = line.rstrip("\n")
                if line.startswith("TITLE="):
                    title = line[6:].strip()
                elif line == "END IONS":
                    if title:
                        stem, scan = split_title(title)
                        key = (stem, scan)
                        if key in wanted:
                            out[key] = parse_peaks(buf)
                    title, buf = None, []
                elif line and line[0].isdigit():
                    buf.append(line)
    return out


def split_title(t):
    s = t.strip()
    if s.lower().endswith(".dta"):
        s = s[:-4]
    parts = s.split(".")
    if len(parts) >= 5 and parts[-1] == "0":
        s = ".".join(parts[:-1])
    parts = s.rsplit(".", 3)
    try:
        return parts[0], int(parts[1])
    except (IndexError, ValueError):
        return s, 0


def parse_peaks(lines):
    peaks = []
    for ln in lines:
        f = ln.split()
        if len(f) >= 2:
            try:
                peaks.append((float(f[0]), float(f[1])))
            except ValueError:
                continue
    return peaks


def coverage(seq, z, peaks):
    """b/y 覆盖度 = 命中离子数 / 理论离子数（±20 ppm）。"""
    ions = b_y_ions(seq, z)
    if not ions or not peaks:
        return 0.0, 0, len(ions)
    mz = sorted(p[0] for p in peaks)
    hit = 0
    j = 0
    import bisect
    for m, _lab in ions:
        tol = m * TOL * 1e-6
        lo, hi = m - tol, m + tol
        k = bisect.bisect_left(mz, lo)
        if k < len(mz) and mz[k] <= hi:
            hit += 1
    return (hit / len(ions) if ions else 0.0), hit, len(ions)


# ---------------------------------------------------------------- 主流程
def main(argv=None):
    ap = argparse.ArgumentParser(description="TP 结果 P13 四件套")
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--skip-peaks", action="store_true", help="跳过 b/y 覆盖度计算")
    args = ap.parse_args(argv)

    root = os.path.realpath(args.workdir)
    out = os.path.realpath(args.out)
    os.makedirs(out, exist_ok=True)

    # ---------- 1) 新肽（P7 表）+ 谱数
    v = pd.read_csv(os.path.join(root, "06b_validate_research", "novel_validated.tsv"),
                    sep="\t", dtype=str, encoding="utf-8-sig")
    v["pep"] = v["peptide"].map(iln)
    keep_peps = set(v.loc[v["final_keep"].astype(str).str.startswith("KEEP"), "pep"])
    # 类别在**肽**级判定：任一行为 KEEP 即保留
    pe_cat = {p: ("保留新肽段" if p in keep_peps else "排除新肽段")
              for p in v["pep"].unique()}

    # 谱数：直接从 E2/P6b 的 PSM 表数（按肽×酶）
    psm = defaultdict(lambda: defaultdict(int))     # pep -> enz -> n
    for sub in ("04_search_disc", "05b_research"):
        base = os.path.join(root, sub, "TP")
        if not os.path.isdir(base):
            continue
        for t in sorted(os.listdir(base)):
            f = os.path.join(base, t, "pFind-Filtered.spectra")
            if not os.path.exists(f):
                continue
            enz = "LL" if "_LL_" in t else "TL"
            for r in load_spectra(f):
                p = iln(r.get("Sequence") or "")
                if p in pe_cat:
                    psm[p][enz] += 1

    # ORF（蛋白）归属
    prot_of = defaultdict(set)
    for p, pr in zip(v["pep"], v["protein"]):
        for tok in str(pr).split("/"):
            tok = tok.strip()
            if tok:
                prot_of[p].add(tok)

    # ---------- 2) 已注释（E1）
    ann = dict(psm=0, pep=set(), prot=set())
    ann_enz = defaultdict(lambda: dict(psm=0, pep=set(), prot=set()))
    abase = os.path.join(root, "04_search_annot", "TP")
    for t in sorted(os.listdir(abase)) if os.path.isdir(abase) else []:
        f = os.path.join(abase, t, "pFind-Filtered.spectra")
        if not os.path.exists(f):
            continue
        enz = "LL" if "_LL_" in t else "TL"
        for r in load_spectra(f):
            prots = [x.strip().rstrip("/") for x in (r.get("Proteins") or "").split(";")]
            l0 = [x for x in prots if x.startswith("L0_")]
            if not l0:
                continue
            p = iln(r.get("Sequence") or "")
            for tgt in (ann, ann_enz[enz]):
                tgt["psm"] += 1
                if p:
                    tgt["pep"].add(p)
                tgt["prot"].update(l0)

    # ---------- 3) 分类统计表
    rows = []
    for enz in ("TL", "LL", "ALL"):
        def n_psm(pep_set):
            if enz == "ALL":
                return sum(sum(psm[p].values()) for p in pep_set)
            return sum(psm[p].get(enz, 0) for p in pep_set)

        if enz == "ALL":
            a_psm, a_pep, a_orf = ann["psm"], len(ann["pep"]), len(ann["prot"])
        else:
            a = ann_enz[enz]
            a_psm, a_pep, a_orf = a["psm"], len(a["pep"]), len(a["prot"])
        for cat in ("保留新肽段", "排除新肽段"):
            # 分酶切时只算在该酶下**有谱**的肽（否则两个酶的数字会一模一样）
            if enz == "ALL":
                ps = [p for p, c in pe_cat.items() if c == cat]
            else:
                ps = [p for p, c in pe_cat.items() if c == cat and psm[p].get(enz, 0) > 0]
            orfs = set()
            for p in ps:
                orfs |= prot_of[p]
            rows.append(dict(Enzyme=ENZ_EN.get(enz, "Combined"), Category=CAT_EN[cat],
                             ORFs=len(orfs), Peptides=len(ps), PSMs=n_psm(ps)))
        rows.append(dict(Enzyme=ENZ_EN.get(enz, "Combined"), Category=CAT_EN["已注释"],
                         ORFs=a_orf, Peptides=a_pep, PSMs=a_psm))
    t1 = pd.DataFrame(rows)
    order = {CAT_EN["已注释"]: 0, CAT_EN["排除新肽段"]: 1, CAT_EN["保留新肽段"]: 2}
    t1["_o"] = t1["Category"].map(order)
    t1 = t1.sort_values(["_o"]).drop(columns=["_o"])
    t1.to_csv(os.path.join(out, "P13_1_category_stats.tsv"), sep="\t",
              index=False, encoding="utf-8-sig")

    # ---------- 4) 质控：相似性 / 覆盖度 / Raw score
    # 每条肽取最优 PSM（Final_Score 最小）
    best = {}
    for sub in ("04_search_disc", "05b_research"):
        base = os.path.join(root, sub, "TP")
        if not os.path.isdir(base):
            continue
        for t in sorted(os.listdir(base)):
            f = os.path.join(base, t, "pFind-Filtered.spectra")
            if not os.path.exists(f):
                continue
            enz = "LL" if "_LL_" in t else "TL"
            for r in load_spectra(f):
                p = iln(r.get("Sequence") or "")
                if p not in pe_cat:
                    continue
                try:
                    fs = float(r.get("Final_Score") or 9e9)
                except ValueError:
                    fs = 9e9
                cur = best.get(p)
                if cur is None or fs < cur["fs"]:
                    best[p] = dict(fs=fs, enz=enz, seq=r.get("Sequence"),
                                   title=r.get("File_Name", ""),
                                   raw=r.get("Raw_Score", ""),
                                   z=r.get("Charge", ""),
                                   mod=r.get("Modification", ""))

    # 峰缓存（只为 QC 覆盖度用）
    peaks = {}
    if not args.skip_peaks:
        wanted = set()
        for p, b in best.items():
            stem, scan = split_title(b["title"])
            wanted.add((stem, scan))
        mgfs = []
        mdir = os.path.join(root, "03_mgf")
        if os.path.isdir(mdir):
            mgfs = [os.path.join(mdir, f) for f in sorted(os.listdir(mdir))
                    if f.endswith(".mgf")]
        print("  读 MGF 取峰：%d 张谱 / %d 个文件" % (len(wanted), len(mgfs)))
        peaks = read_mgf_peaks(mgfs, wanted)
        print("  命中峰谱 %d / %d" % (len(peaks), len(wanted)))

    qc = []
    cos_of = {}
    for _, r in v.iterrows():
        if r["v_best_cos"] not in (None, ""):
            try:
                c = float(r["v_best_cos"])
            except ValueError:
                continue
            p = r["pep"]
            if p not in cos_of or c > cos_of[p]:
                cos_of[p] = c
    for p, cat in pe_cat.items():
        b = best.get(p, {})
        seq = b.get("seq") or p
        try:
            z = int(float(b.get("z") or 2))
        except ValueError:
            z = 2
        stem, scan = split_title(b.get("title", ""))
        pk = peaks.get((stem, scan), [])
        cov, hit, tot = coverage(seq, z, pk) if pk else (float("nan"), 0, 0)
        try:
            raw = float(b.get("raw") or "nan")
        except ValueError:
            raw = float("nan")
        qc.append(dict(peptide=p, category=CAT_EN[cat],
                       enzyme=ENZ_NAME.get(b.get("enz", ""), ""),
                       cosine=cos_of.get(p, float("nan")),
                       ion_coverage=cov, ion_hit=hit, ion_total=tot,
                       raw_score=raw, charge=z,
                       n_spectra=sum(psm[p].values()),
                       file=b.get("title", "")))
    qc = pd.DataFrame(qc)
    qc.to_csv(os.path.join(out, "P13_2_qc_per_peptide.tsv"), sep="\t",
              index=False, encoding="utf-8-sig")

    print("\n== 分类统计 ==")
    print(t1.to_string(index=False))
    print("\n== 质控（中位数）==")
    g = qc.groupby(["category", "enzyme"]).agg(
        n=("peptide", "size"), cosine=("cosine", "median"),
        ion_coverage=("ion_coverage", "median"), raw_score=("raw_score", "median"),
        spectra=("n_spectra", "median"))
    print(g.to_string())
    return 0


if __name__ == "__main__":
    sys.exit(main())
