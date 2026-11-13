# -*- coding: utf-8 -*-
"""TP P13 汇总表（v2，口径统一到"唯一位点"）。

为什么重算
----------
旧表把 ORF 数按 `protein` 字段直接并集，而 139 条保留肽各自同时记在
两个 ORF 名下（L1 ab initio 与 L2 六框对同一物理位点各给一个名字）。
逐条比对坐标：139 对全部同链且重叠 >80%，即同一个位点 → 旧值 114/115 是重复计数。

本脚本统一口径
--------------
  GSSP / 肽段 = 唯一序列
  ORF         = 唯一**位点**，用事件映射给出的主 ORF（event_orf）代表；
                同一肽的第二个同座位名已被去重
  酶切        = 该肽有 PSM 的任务（TL→Trypsin，LL→LysargiNase），一条肽只归一酶

输出
----
  P13_6_table_v2.tsv    模板表（Type / Parameters / Trypsin / LysargiNase / Total）
  P13_6_tiers_v2.tsv    证据分级明细
  表本身同时打印到 stdout（markdown），由调用方落盘
"""
import argparse
import os
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))

import pandas as pd

from orfscape.core.novel import iln, load_spectra

ENZ = {"TL": "Trypsin", "LL": "LysargiNase"}
SUB_CN = {"intergenic_only": "Intergenic spacer (ITS)",
          "same_strand_other_frame": "Internal alternative frame",
          "antisense_overlap": "Antisense translation",
          "nterm_extension": "N-termini correction"}
NOVEL_SUBS = ["intergenic_only", "same_strand_other_frame", "antisense_overlap"]


def core(p):
    """镜像核心：去两端 K/R。"""
    s = str(p)
    while s and s[0] in "KR":
        s = s[1:]
    while s and s[-1] in "KR":
        s = s[:-1]
    return s


def _inside(base, path):
    """path 必须落在 base 之内（两者均已 realpath 规范化）。"""
    return path == base or path.startswith(base + os.sep)


def spectra_counts(root, pep_set):
    """肽 × 酶 → PSM 数（E2 + P6b）。"""
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
                if p in pep_set:
                    psm[p][enz] += 1
    return psm


def annotated_counts(root):
    """E1 里命中 L0_* 的 PSM / 肽 / ORF，分酶 + 并集。"""
    tot = dict(psm=0, pep=set(), orf=set())
    per = {e: dict(psm=0, pep=set(), orf=set()) for e in ("TL", "LL")}
    base = os.path.join(root, "04_search_annot", "TP")
    for t in sorted(os.listdir(base)) if os.path.isdir(base) else []:
        f = os.path.join(base, t, "pFind-Filtered.spectra")
        if not os.path.exists(f):
            continue
        enz = "LL" if "_LL_" in t else "TL"
        for r in load_spectra(f):
            prots = [x.strip().rstrip("/") for x in (r.get("Proteins") or "").split(";")]
            l0 = [x for x in prots if x.startswith("L0_")]
            if not l0:
                continue
            p = iln(r.get("Sequence") or "")
            for tgt in (tot, per[enz]):
                tgt["psm"] += 1
                if p:
                    tgt["pep"].add(p)
                tgt["orf"].update(l0)
    return tot, per


def _overlap_same(a, b, frac=0.8):
    """两个 ORF 区间是否同 contig、同链、重叠 ≥ frac（取较短者作分母）。"""
    if not a or not b:
        return False
    if a[0] != b[0] or a[1] != b[1]:
        return False
    lo, hi = max(a[2], b[2]), min(a[3], b[3])
    if hi < lo:
        return False
    shorter = min(a[3] - a[2] + 1, b[3] - b[2] + 1)
    return (hi - lo + 1) / max(1, shorter) >= frac


def locus_components(om, prot_sets):
    """把"同一物理位点、不同预测名"的 ORF 合并。

    prot_sets: {肽: {ORF, ...}}，同一肽的多个候选 ORF 若同座位则连边；
    最后取连通分量数 = 唯一位点数。
    """
    parent = {}

    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    span = {}
    for orfs in prot_sets.values():
        for o in orfs:
            if o not in span:
                span[o] = om.orf_span(o)
            find(o)
        lst = sorted(orfs)
        for i in range(len(lst)):
            for j in range(i + 1, len(lst)):
                if _overlap_same(span[lst[i]], span[lst[j]]):
                    union(lst[i], lst[j])
    return {o: find(o) for o in parent}


def main(argv=None):
    ap = argparse.ArgumentParser(description="TP P13 汇总表 v2")
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--res", required=True, help="P13_results 目录")
    args = ap.parse_args(argv)

    root = os.path.realpath(args.workdir)
    res = os.path.realpath(args.res)
    if not _inside(root, res):
        raise SystemExit("输出目录必须位于项目目录内：%s" % res)

    v = pd.read_csv(os.path.join(root, "06b_validate_research", "novel_validated.tsv"),
                    sep="\t", dtype=str, encoding="utf-8-sig")
    v["pep"] = v["peptide"].map(iln)
    keep = v[v["final_keep"].astype(str).str.startswith("KEEP")]
    peps = keep.drop_duplicates("pep")[["pep", "peptide", "enzymes_str"]].copy()
    peps = peps.rename(columns={"enzymes_str": "enz"})
    assert peps["pep"].is_unique

    ev = pd.read_csv(os.path.join(res, "P13_3_events_per_peptide.tsv"),
                     sep="\t", dtype=str, encoding="utf-8-sig")
    ev["pep"] = ev["peptide"].map(iln)
    e1 = ev.drop_duplicates("pep").set_index("pep")
    peps["orf"] = peps["pep"].map(lambda p: e1.loc[p, "event_orf"])
    # 注意：事件表里 `sub` 是中文标签，`event_sub` 才是英文键
    peps["sub"] = peps["pep"].map(lambda p: e1.loc[p, "event_sub"])

    loci = keep["protein"].nunique()
    peps_multi = int((keep.groupby("pep")["protein"].nunique() > 1).sum())

    # 谱数：扫全部新肽（保留 + 排除都要），否则排除组的谱数会是 0
    psm = spectra_counts(root, set(v["pep"]))
    ann_tot, ann_per = annotated_counts(root)

    rows = []

    def add(sec, par, t, l, tot):
        rows.append(dict(Section=sec, Parameters=par,
                         Trypsin=int(t), LysargiNase=int(l), Total=int(tot)))

    add("Total identification", "PSMs",
        ann_per["TL"]["psm"], ann_per["LL"]["psm"], ann_tot["psm"])
    add("Total identification", "Peptides",
        len(ann_per["TL"]["pep"]), len(ann_per["LL"]["pep"]), len(ann_tot["pep"]))
    add("Total identification", "ORFs",
        len(ann_per["TL"]["orf"]), len(ann_per["LL"]["orf"]), len(ann_tot["orf"]))

    def n_psm(pep_list, enz=None):
        if enz:
            return sum(psm[p].get(enz, 0) for p in pep_list)
        return sum(sum(psm[p].values()) for p in pep_list)

    t_pep = peps[peps["enz"] == "TL"]
    l_pep = peps[peps["enz"] == "LL"]
    add("After filtering", "PSMs",
        n_psm(list(t_pep["pep"]), "TL"), n_psm(list(l_pep["pep"]), "LL"),
        n_psm(list(peps["pep"])))
    add("After filtering", "Peptides", len(t_pep), len(l_pep), len(peps))
    add("After filtering", "ORFs",
        t_pep["orf"].nunique(), l_pep["orf"].nunique(), peps["orf"].nunique())

    for sub in ["nterm_extension"] + NOVEL_SUBS:
        s = peps[peps["sub"] == sub]
        p_t = int((s["enz"] == "TL").sum())
        p_l = int((s["enz"] == "LL").sum())
        add("Reannotation events", SUB_CN[sub] + " | GSSPs", p_t, p_l, p_t + p_l)
        add("Reannotation events", SUB_CN[sub] + " | ORFs",
            s[s["enz"] == "TL"]["orf"].nunique(),
            s[s["enz"] == "LL"]["orf"].nunique(), s["orf"].nunique())
    nov = peps[peps["sub"].isin(NOVEL_SUBS)]
    add("Reannotation events", "Novel ORF (union) | GSSPs",
        int((nov["enz"] == "TL").sum()), int((nov["enz"] == "LL").sum()), len(nov))
    add("Reannotation events", "Novel ORF (union) | ORFs",
        nov[nov["enz"] == "TL"]["orf"].nunique(),
        nov[nov["enz"] == "LL"]["orf"].nunique(), nov["orf"].nunique())
    add("Reannotation events", "All events (union) | ORFs", 0, 0, peps["orf"].nunique())

    # 证据分级
    grp = defaultdict(lambda: dict(pep=set(), t=set(), l=set()))
    for _, r in peps.iterrows():
        g = grp[r["orf"]]
        g["pep"].add(r["pep"])
        g["t" if r["enz"] == "TL" else "l"].add(r["pep"])

    def tier(g):
        if {core(p) for p in g["t"]} & {core(p) for p in g["l"]}:
            return "1"
        if g["t"] and g["l"]:
            return "2"
        if g["t"]:
            return "3.1-multi" if len(g["t"]) >= 2 else "3.1-single"
        if g["l"]:
            return "3.2-multi" if len(g["l"]) >= 2 else "3.2-single"
        return "?"

    cnt = defaultdict(lambda: [0, 0, 0])
    for orf, g in grp.items():
        c = cnt[tier(g)]
        c[0] += 1
        c[1] += len(g["pep"])
        c[2] += sum(sum(psm[p].values()) for p in g["pep"])

    tier_rows = []
    for k, name in [("1", "Supported by at least one mirror GSSP pair"),
                    ("2", "Supported by two proteases without mirror GSSP"),
                    ("3.1-multi", "exclusive Trypsin, multi-peptide"),
                    ("3.1-single", "exclusive Trypsin, single-peptide"),
                    ("3.2-multi", "exclusive LysargiNase, multi-peptide"),
                    ("3.2-single", "exclusive LysargiNase, single-peptide")]:
        tier_rows.append(dict(Type="Novel ORF", Parameters=name, ORF=cnt[k][0],
                              Peptides=cnt[k][1], PSMs=cnt[k][2]))

    # ---------- ORF 互斥归类（用户 2026-11-11 定）
    # 同一个 ORF 可能承载多个亚型的肽；为了让 ORF 列可相加，每个 ORF 只归一类：
    #   N 端校正 > 不同框 > 反义链 > 基因间区
    # 理由：N 端校正是对已知基因的定点修正（最强主张）；不同框/反义链都表示与已知基因重叠，
    # 基因间区是"无重叠"兜底，一旦有重叠证据就不再算基因间区。
    PRIO = {"nterm_extension": 0, "same_strand_other_frame": 1,
            "antisense_overlap": 2, "intergenic_only": 3}
    orf_sub = defaultdict(set)
    for _, r in peps.iterrows():
        orf_sub[r["orf"]].add(r["sub"])
    orf_cat = {o: min(s, key=lambda x: PRIO[x]) for o, s in orf_sub.items()}
    peps["cat"] = peps["orf"].map(orf_cat)
    moved = {o: sorted(s) for o, s in orf_sub.items() if len(s) > 1}

    v3 = []
    for sub in ["nterm_extension"] + NOVEL_SUBS:
        s = peps[peps["cat"] == sub]
        v3.append(dict(Category=SUB_CN[sub], ORF=s["orf"].nunique(),
                       ORF_T=s[s["enz"] == "TL"]["orf"].nunique(),
                       ORF_L=s[s["enz"] == "LL"]["orf"].nunique(),
                       GSSP_by_ORF=len(s),
                       GSSP_T_by_ORF=int((s["enz"] == "TL").sum()),
                       GSSP_L_by_ORF=int((s["enz"] == "LL").sum()),
                       GSSP_by_peptide=int((peps["sub"] == sub).sum()),
                       GSSP_T_by_peptide=int(((peps["sub"] == sub) &
                                              (peps["enz"] == "TL")).sum()),
                       GSSP_L_by_peptide=int(((peps["sub"] == sub) &
                                              (peps["enz"] == "LL")).sum()),
                       PSMs_by_ORF=n_psm(list(s["pep"])),
                       PSMs_T_by_ORF=n_psm(list(s[s["enz"] == "TL"]["pep"]), "TL"),
                       PSMs_L_by_ORF=n_psm(list(s[s["enz"] == "LL"]["pep"]), "LL"),
                       PSMs_by_peptide=n_psm(list(peps[peps["sub"] == sub]["pep"]))))
    v3 = pd.DataFrame(v3)
    tot = dict(Category="Total", ORF=int(v3["ORF"].sum()),
               ORF_T=int(v3["ORF_T"].sum()), ORF_L=int(v3["ORF_L"].sum()),
               PSMs_by_ORF=int(v3["PSMs_by_ORF"].sum()),
               PSMs_T_by_ORF=int(v3["PSMs_T_by_ORF"].sum()),
               PSMs_L_by_ORF=int(v3["PSMs_L_by_ORF"].sum()),
               PSMs_by_peptide=int(v3["PSMs_by_peptide"].sum()))
    for c in ["GSSP_by_ORF", "GSSP_T_by_ORF", "GSSP_L_by_ORF",
              "GSSP_by_peptide", "GSSP_T_by_peptide", "GSSP_L_by_peptide"]:
        tot[c] = int(v3[c].sum())
    v3.loc[len(v3)] = tot
    v3.to_csv(os.path.join(res, "P13_6_table_v3.tsv"), sep="\t", index=False,
              encoding="utf-8-sig")

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(res, "P13_6_table_v2.tsv"), sep="\t", index=False,
              encoding="utf-8-sig")
    pd.DataFrame(tier_rows).to_csv(os.path.join(res, "P13_6_tiers_v2.tsv"), sep="\t",
                                   index=False, encoding="utf-8-sig")

    def gv(sec, par, col):
        r = df[(df.Section == sec) & (df.Parameters == par)]
        return int(r.iloc[0][col]) if len(r) else 0

    L = ["| Type | Parameters | Trypsin | LysargiNase | Total |",
         "|---|---|---|---|---|"]
    for sec in ["Total identification", "After filtering"]:
        for par in ["PSMs", "Peptides", "ORFs"]:
            r = df[(df.Section == sec) & (df.Parameters == par)].iloc[0]
            L.append("| %s | %s | %s | %s | %s |" % (
                sec if par == "PSMs" else "", par,
                format(int(r.Trypsin), ","), format(int(r.LysargiNase), ","),
                format(int(r.Total), ",")))
    L.append("|  |  |  |  |  |")
    for sub in ["nterm_extension"] + NOVEL_SUBS:
        sec = "Reannotation events"
        pt = gv(sec, SUB_CN[sub] + " | GSSPs", "Trypsin")
        pl = gv(sec, SUB_CN[sub] + " | GSSPs", "LysargiNase")
        oa = gv(sec, SUB_CN[sub] + " | ORFs", "Total")
        nm = "N-termini correction" if sub == "nterm_extension" else SUB_CN[sub]
        L.append("| %s | %s | %d (%d/%d) | | %d |" % (
            "Reannotation events" if sub == "nterm_extension" else "", nm,
            pt + pl, pt, pl, oa))
    L.append("|  |  |  |  |  |")
    L.append("| Novel ORF | Supported by at least one mirror GSSP pair | | | %d |"
             % cnt["1"][0])
    L.append("| | Supported by two proteases without mirror GSSP | | | %d |" % cnt["2"][0])
    L.append("| | Supported exclusively by Trypsin GSSPs | %d | | %d |"
             % (cnt["3.1-multi"][0] + cnt["3.1-single"][0],
                cnt["3.1-multi"][0] + cnt["3.1-single"][0]))
    L.append("| | Supported exclusively by LysargiNase GSSPs | | %d | %d |"
             % (cnt["3.2-multi"][0] + cnt["3.2-single"][0],
                cnt["3.2-multi"][0] + cnt["3.2-single"][0]))
    print("\n".join(L))
    print("\n== 证据分级 ==")
    print(pd.DataFrame(tier_rows).to_string(index=False))
    print("\n旧口径 ORF（protein 名并集，含同座位重名）: %d ｜ 同肽双 ORF 的肽数: %d"
          % (loci, peps_multi))
    print("事件 ORF 并集: %d" % peps["orf"].nunique())

    # ---------- 四、三类别的"唯一位点"ORF 数（含排除组，用于汇总表自洽）
    from tp_p13_pdf_inputs import OrfMap
    om = OrfMap(root)

    def prot_sets(df):
        out = defaultdict(set)
        for p, pr in zip(df["pep"], df["protein"]):
            for tok in str(pr).split("/"):
                tok = tok.strip()
                if tok:
                    out[p].add(tok)
        return out

    keep_ps = prot_sets(keep)
    excl_ps = prot_sets(v[~v["pep"].isin(set(keep["pep"]))])
    comp_k = locus_components(om, keep_ps)
    comp_e = locus_components(om, excl_ps)

    def n_locus(comp, ps):
        return len({comp[o] for orfs in ps.values() for o in orfs if o in comp})

    out = pd.DataFrame([
        dict(Category="Annotated", ORF_names=len(ann_tot["orf"]),
             Unique_loci=len(ann_tot["orf"]), Peptides=len(ann_tot["pep"]), PSMs=ann_tot["psm"]),
        dict(Category="Filtered", ORF_names=len({o for s in excl_ps.values() for o in s}),
             Unique_loci=n_locus(comp_e, excl_ps), Peptides=len(excl_ps),
             PSMs=n_psm(list(excl_ps))),
        dict(Category="Obtained", ORF_names=len({o for s in keep_ps.values() for o in s}),
             Unique_loci=n_locus(comp_k, keep_ps), Peptides=len(keep_ps),
             PSMs=n_psm(list(keep_ps))),
    ])
    out.to_csv(os.path.join(res, "P13_6_category_v2.tsv"), sep="\t", index=False,
               encoding="utf-8-sig")
    print("\n== 三类别的 ORF 口径对照 ==")
    print(out.to_string(index=False))

    print("\n== ORF 互斥归类（N 端校正 > 不同框 > 反义链 > 基因间区）==")
    print(pd.DataFrame(v3).to_string(index=False))
    print("\n承载多个亚型的 ORF：")
    for o, s in moved.items():
        print("  %-32s %s → 归 %s（%d 肽）"
              % (o, "/".join(s), orf_cat[o], int((peps["orf"] == o).sum())))
    return 0


if __name__ == "__main__":
    sys.exit(main())
