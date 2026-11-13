# -*- coding: utf-8 -*-
"""Novel ORF 的证据分级（用户 2026-11-11 口径）。

   1 级  有镜像肽段证据：该 ORF 内存在"跨酶同核心肽"（去两端 K/R 后相同，T 与 L 各见一次）
   2 级  无镜像，但有**两种酶切**证据：该 ORF 的肽同时来自 Trypsin 与 LysargiNase
   3.1  只有 Trypsin 肽：再分 多肽段（≥2）/ 单肽段（=1）
   3.2  只有 LysargiNase 肽：同上

输出
  P13_5_orf_tiers.tsv          每个 ORF 一行（级别、肽/谱数、酶、镜像核心）
  P13_5_orf_tier_summary.tsv   级别 × ORF/肽段/谱图
"""
import argparse
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pandas as pd

from orfscape.core.novel import iln, load_spectra

ENZ = {"TL": "Trypsin", "LL": "LysargiNase"}


def core(p):
    """镜像核心：去两端 K/R。"""
    s = str(p)
    while s and s[0] in "KR":
        s = s[1:]
    while s and s[-1] in "KR":
        s = s[:-1]
    return s


def main(argv=None):
    ap = argparse.ArgumentParser(description="Novel ORF 证据分级")
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--res", required=True, help="P13_results 目录")
    args = ap.parse_args(argv)

    root = os.path.realpath(args.workdir)
    res = os.path.realpath(args.res)

    v = pd.read_csv(os.path.join(root, "06b_validate_research", "novel_validated.tsv"),
                    sep="\t", dtype=str, encoding="utf-8-sig")
    v["pep"] = v["peptide"].map(iln)
    keep = v[v["final_keep"].astype(str).str.startswith("KEEP")].copy()

    ev = pd.read_csv(os.path.join(res, "P13_3_events_per_peptide.tsv"),
                     sep="\t", dtype=str, encoding="utf-8-sig")
    ev["pep"] = ev["peptide"].map(iln)
    orf_of = dict(zip(ev["pep"], ev["event_orf"].fillna("")))

    # 每肽谱数（分酶），直接数 PSM
    psm = defaultdict(lambda: defaultdict(int))
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
                if p:
                    psm[p][enz] += 1

    # ---- 按 ORF 归并
    by_orf = defaultdict(lambda: dict(peps={}, enz=defaultdict(set), n=0))
    for _, r in keep.iterrows():
        p = r["pep"]
        orf = orf_of.get(p) or str(r["protein"]).split("/")[0]
        rec = by_orf[orf]
        rec["peps"][p] = r["peptide"]
        # novel_validated 里的酶代号是任务名（TL=Trypsin / LL=LysargiNase），
        # 也可能写成 "TL/LL"（两个酶都测到）
        for e in str(r.get("enzymes_str", "")).split("/"):
            e = e.strip().upper()
            if e == "TL":
                rec["enz"]["T"].add(p)
            elif e == "LL":
                rec["enz"]["L"].add(p)

    rows = []
    for orf, rec in sorted(by_orf.items()):
        peps = set(rec["peps"])
        t_set, l_set = rec["enz"]["T"], rec["enz"]["L"]
        core_t = {core(p) for p in t_set}
        core_l = {core(p) for p in l_set}
        mirror = core_t & core_l
        if mirror:
            lvl = "1 Mirror GSSP pair"
        elif t_set and l_set:
            lvl = "2 Two proteases (no mirror)"
        elif t_set:
            lvl = "3.1 Trypsin only, multi-peptide" if len(t_set) >= 2 else "3.1 Trypsin only, single-peptide"
        elif l_set:
            lvl = "3.2 LysargiNase only, multi-peptide" if len(l_set) >= 2 else "3.2 LysargiNase only, single-peptide"
        else:
            lvl = "Unassigned"
        n_psm = sum(sum(psm[p].values()) for p in peps)
        rows.append(dict(ORF=orf, level=lvl, Peptides=len(peps), PSMs=int(n_psm),
                         Trypsin_peptides=len(t_set), LysargiNase_peptides=len(l_set),
                         mirror_core=";".join(sorted(mirror)) if mirror else "",
                         peptide_list=";".join(sorted(peps))))
    df = pd.DataFrame(rows)
    order = ["1 Mirror GSSP pair", "2 Two proteases (no mirror)",
             "3.1 Trypsin only, multi-peptide", "3.1 Trypsin only, single-peptide",
             "3.2 LysargiNase only, multi-peptide", "3.2 LysargiNase only, single-peptide", "Unassigned"]
    df["_o"] = df["level"].map({k: i for i, k in enumerate(order)})
    df = df.sort_values(["_o", "ORF"]).drop(columns=["_o"])
    df.to_csv(os.path.join(res, "P13_5_orf_tiers.tsv"), sep="\t", index=False,
              encoding="utf-8-sig")

    summ = (df.groupby("level").agg(ORFs=("ORF", "size"), Peptides=("Peptides", "sum"),
                               PSMs=("PSMs", "sum")).reset_index())
    summ["_o"] = summ["level"].map({k: i for i, k in enumerate(order)})
    summ = summ.sort_values("_o").drop(columns=["_o"])
    summ.to_csv(os.path.join(res, "P13_5_orf_tier_summary.tsv"), sep="\t",
                index=False, encoding="utf-8-sig")
    print(summ.to_string(index=False))
    print("\nORFs %d | peptides %d" % (len(df), int(df["Peptides"].sum())))
    return 0


if __name__ == "__main__":
    sys.exit(main())
