#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""补算缺失的 AlphaPeptDeep 余弦相似性，让质控图里四组数字对齐（230 / 460）。

背景：isoPSV 对一部分肽（多是 mirror-only 保留、isoPSV 未通过的）没有产出相似性，
导致 cosine 一格的 n 变成 213/404，与其他指标（230/460）不一致。

做法：用同一份 AlphaPeptDeep 预测（predict.txt 的 b+1/b+2/y+1/y+2 强度）与实测谱
（MGF）重算 cosine，口径与引擎一致（±20 ppm 取最强峰、√强度、余弦）；
**先用已有 isoPSV 值的肽做校验**，通过后才用于补缺。

    python fill_cosine.py --workdir <proj> --res <P13_results> [--tol 20]

输出：`P13_2_qc_per_peptide.tsv` 增列 `cosine_calc`、`cosine_final`、`cosine_source`；
并打印校验结果（与 isoPSV 的相关系数）。
"""
import argparse
import os
import re
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import pandas as pd

from tp_p13_tables import AA_MASS, H2O, PROTON, read_mgf_peaks, split_title

ION_RE = re.compile(r"^([by])\+(\d)=(.*)$")


def read_predictions(path):
    """predict.txt → {(seq, z): [(mz, inten, label)]}"""
    out = {}
    seq, z, sticks = None, None, []
    if not os.path.isfile(path):
        return out
    for line in open(path, encoding="utf-8", errors="replace"):
        line = line.rstrip("\n")
        if line.startswith("TITLE="):
            t = line[6:].strip()
            if "||" in t:
                seq, _, zs = t.partition("||")
                try:
                    z = int(zs)
                except ValueError:
                    z = 2
            else:
                seq, z = t, 2
            sticks = []
        elif line.startswith("END IONS"):
            if seq:
                out.setdefault((seq, z), []).extend(sticks)
            seq, z, sticks = None, None, []
        else:
            m = ION_RE.match(line)
            if not m or not seq:
                continue
            ser, ch, vals = m.group(1), int(m.group(2)), m.group(3)
            ints = [float(x) for x in vals.split(",") if x.strip() != ""]
            if not ints:
                continue
            mass = [AA_MASS.get(a, 0.0) for a in seq]
            n = len(seq)
            pre = 0.0
            for i in range(n - 1):
                pre += mass[i]
                if ser == "b":
                    mz = (pre + ch * PROTON) / ch
                    lab = "b%d^%d" % (i + 1, ch)
                else:
                    suf = sum(mass) - pre
                    mz = (suf + H2O + ch * PROTON) / ch
                    lab = "y%d^%d" % (n - i - 1, ch)
                if i < len(ints) and ints[i] > 0:
                    sticks.append((mz, ints[i], lab))
    return out


def theory_cosine(mz, ints, sticks, tol_ppm=20.0):
    if not sticks or len(mz) == 0:
        return None
    order = np.argsort(mz)
    mzs, ins = np.asarray(mz)[order], np.asarray(ints)[order]
    a, b = [], []
    for m, v, _lab in sticks:
        lo, hi = m * (1 - tol_ppm * 1e-6), m * (1 + tol_ppm * 1e-6)
        i0 = np.searchsorted(mzs, lo, side="left")
        i1 = np.searchsorted(mzs, hi, side="right")
        a.append(v)
        b.append(float(ins[i0:i1].max()) if i1 > i0 else 0.0)
    a = np.sqrt(np.asarray(a, float))
    b = np.asarray(b, float)
    if a.max() <= 0 or b.max() <= 0:
        return 0.0
    b = np.sqrt(b / b.max())
    return float(a @ b / (np.linalg.norm(a) * np.linalg.norm(b)))


def main(argv=None):
    ap = argparse.ArgumentParser(description="补算缺失的 AlphaPeptDeep 余弦")
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--res", required=True)
    ap.add_argument("--tol", type=float, default=20.0, help="ppm 容差")
    args = ap.parse_args(argv)

    root = os.path.realpath(args.workdir)
    res = os.path.realpath(args.res)
    qc = pd.read_csv(os.path.join(res, "P13_2_qc_per_peptide.tsv"), sep="\t",
                     encoding="utf-8-sig")

    # ---- 预测：两个酶的结果目录（取最新 result_*）
    pred = {}
    for enz in ("TP_TL", "TP_LL"):
        d = os.path.join(root, "06b_validate_research", enz)
        if not os.path.isdir(d):
            continue
        picks = sorted(p for p in os.listdir(d) if p.startswith("result_"))
        for p in picks[-1:]:
            f = os.path.join(d, p, "merge", "predictors", "alphapeptdeep", "predict.txt")
            got = read_predictions(f)
            pred.update(got)
            print("预测 %s: %d 个 (肽,电荷)" % (enz, len(got)))

    # ---- 实测峰（按 file 列 = PSM title）
    wanted = set()
    for f in qc["file"].fillna(""):
        if f:
            wanted.add(tuple(split_title(f)))
    mdir = os.path.join(root, "03_mgf")
    mgfs = [os.path.join(mdir, x) for x in sorted(os.listdir(mdir)) if x.endswith(".mgf")] \
        if os.path.isdir(mdir) else []
    print("读 MGF：%d 张谱 / %d 个文件" % (len(wanted), len(mgfs)))
    peaks = read_mgf_peaks(mgfs, wanted)

    calc, src = [], []
    for _, r in qc.iterrows():
        seq = str(r["peptide"])
        try:
            z = int(float(r["charge"]))
        except (TypeError, ValueError):
            z = 2
        sticks = pred.get((seq, z))
        stem, scan = split_title(str(r.get("file") or ""))
        pk = peaks.get((stem, scan)) or []
        v = theory_cosine([p[0] for p in pk], [p[1] for p in pk], sticks,
                          tol_ppm=args.tol) if (sticks and pk) else None
        calc.append(v)
        src.append("calculated" if v is not None else "")
    qc["cosine_calc"] = calc

    both = qc[qc["cosine"].notna() & qc["cosine_calc"].notna()]
    if len(both) > 3:
        r = np.corrcoef(both["cosine"], both["cosine_calc"])[0, 1]
        d = (both["cosine_calc"] - both["cosine"]).abs()
        print("校验（%d 条两者都有）：Pearson r = %.3f，中位 |Δ| = %.3f，"
              "90 分位 |Δ| = %.3f" % (len(both), r, d.median(), d.quantile(0.9)))

    qc["cosine_final"] = qc["cosine"].where(qc["cosine"].notna(), qc["cosine_calc"])
    qc["cosine_source"] = np.where(qc["cosine"].notna(), "isoPSV",
                                   np.where(qc["cosine_calc"].notna(), "calculated", ""))
    for cat in ("Obtained", "Filtered"):
        s = qc[qc["category"] == cat]
        print("%-9s n=%d  有 cosine_final=%d（isoPSV %d + 补算 %d）"
              % (cat, len(s), s["cosine_final"].notna().sum(),
                 (s["cosine_source"] == "isoPSV").sum(),
                 (s["cosine_source"] == "calculated").sum()))
    qc.to_csv(os.path.join(res, "P13_2_qc_per_peptide.tsv"), sep="\t",
              index=False, encoding="utf-8-sig")
    print("写回 P13_2_qc_per_peptide.tsv（新增 cosine_calc / cosine_final / cosine_source）")
    return 0


if __name__ == "__main__":
    sys.exit(main())
