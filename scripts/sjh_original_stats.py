# -*- coding: utf-8 -*-
"""SJH 原始搜库结果（isoPSV 流程）的肽段层面统计。

数据来源
--------
`original_result/result_Trypsin(2)/result_*/1.unlabeled_result/2.novel_peptide_competition_report.tsv`
`original_result/result_R29/result_*/1.unlabeled_result/2.novel_peptide_competition_report.tsv`
——isoPSV 的"新肽竞争报告"，一行 = 一张谱，`a.Sequence` 是原始搜库（ORF 库）提出的
新肽候选，`credible_novel` / `novel_status` 是 isoPSV 的过滤结论。

三级口径
--------
1. **所有新肽段（识别）**：竞争报告里出现的全部新肽候选（去重、I/L 归一）
2. **过滤新肽段（isoPSV 过滤）**：`credible_novel=TRUE` 的（通过六重存疑规则）
3. **保留新肽段（人工清单）**：SJH 提供的 498 条清单，映射回原始搜库的 1,282 张谱

用法
----
    python scripts/sjh_original_stats.py --original-dir <original_result 目录> \
        --table1 <Table1_novel_peptide_QC.tsv> --out <输出目录>
"""
import argparse
import glob
import os
from collections import Counter

import pandas as pd


def iln(s):
    return str(s or "").strip().upper().replace("I", "L")


def _read_faa(path):
    recs, cur, buf = {}, None, []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if cur is not None:
                    recs[cur] = "".join(buf)
                cur = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if cur is not None:
        recs[cur] = "".join(buf)
    return recs


def find_report(original_dir, which):
    """which: 'Trypsin' / 'R29' → 该流程目录下最新的竞争报告"""
    pats = [os.path.join(original_dir, "result_%s*" % which, "result_*",
                         "1.unlabeled_result",
                         "2.novel_peptide_competition_report.tsv")]
    hits = [p for pat in pats for p in glob.glob(pat)]
    if not hits:
        # 目录名带括号的兼容（result_Trypsin(2)）
        pats = [os.path.join(original_dir, "*%s*" % which, "result_*",
                             "1.unlabeled_result",
                             "2.novel_peptide_competition_report.tsv")]
        hits = [p for pat in pats for p in glob.glob(pat)]
    return sorted(hits)[-1] if hits else ""


def find_spectra(original_dir, which):
    pats = [os.path.join(original_dir, which, "*pFind-Filtered.spectra")]
    hits = [p for pat in pats for p in glob.glob(pat)]
    return sorted(hits)[-1] if hits else ""


def con_map(original_dir):
    """File_Name → (该谱是否只命中污染/自切, 主 ORF 名)。

    竞争报告本身不带蛋白名（`a.Proteins` 全空），但它的行与原始 pFind 表一一对应，
    用 File_Name 就能①剔掉污染/自切的谱、②取回该谱对应的 ORF（计数用）。
    """
    out = {}
    for which in ("Trypsin", "R29"):
        p = find_spectra(original_dir, which)
        if not p:
            continue
        x = pd.read_csv(p, sep="\t", dtype=str, encoding="utf-8-sig",
                        encoding_errors="replace", usecols=["File_Name", "Proteins"])
        for fn, pr in zip(x["File_Name"], x["Proteins"].fillna("")):
            toks = [s.strip() for s in str(pr).split("/")
                    if s.strip().startswith("orf|")]
            out[fn] = (str(pr).startswith("CON_"), toks[0] if toks else "")
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description="isoPSV 流程的肽段层面统计")
    ap.add_argument("--original-dir", required=True)
    ap.add_argument("--table1", required=True, help="498 条清单（Table1_novel_peptide_QC.tsv）")
    ap.add_argument("--annotation", default="",
                    help="注释蛋白 fasta（可选）：用来旁证候选里有序列其实在注释中的")
    ap.add_argument("--out", required=True)
    args = ap.parse_args(argv)

    runs = []
    for which, tag in (("Trypsin", "Trypsin"), ("R29", "R29")):
        rep = find_report(args.original_dir, which)
        if not rep:
            raise SystemExit("找不到 %s 的竞争报告" % which)
        d = pd.read_csv(rep, sep="\t", dtype=str, encoding="utf-8-sig",
                        encoding_errors="replace")
        d["enzyme"] = tag
        d["peptide"] = d["a.Sequence"].fillna("").map(iln)
        d["credible"] = d["credible_novel"].astype(str).str.upper().eq("TRUE")
        runs.append(d)
        print("%-8s 行 %6d | 唯一候选 %5d | credible 行 %6d / 肽 %5d | %s"
              % (tag, len(d), d["peptide"].nunique(), int(d["credible"].sum()),
                 d.loc[d["credible"], "peptide"].nunique(), os.path.basename(rep)))
    d = pd.concat(runs, ignore_index=True)
    d = d[d["peptide"] != ""]
    # 剔掉污染/自切的谱（否则"所有新肽段"里混进上千条污染肽），并取回每张谱的 ORF
    cmap = con_map(args.original_dir)
    d["is_con"] = d["File_Name"].map(lambda f: cmap.get(f, (False, ""))[0])
    d["orf"] = d["File_Name"].map(lambda f: cmap.get(f, (False, ""))[1])
    n_con_rows = int(d["is_con"].sum())
    d = d[~d["is_con"]].copy()

    # ---- 498 清单
    t1 = pd.read_csv(args.table1, sep="\t", dtype=str, encoding="utf-8-sig")
    t1["peptide"] = t1["peptide"].map(iln)
    # Table1 的酶列是 T/L 缩写，这里映射成流程名，便于按酶分组取 ORF
    _E = {"T": "Trypsin", "L": "R29"}
    t1["_enz"] = t1.get("enzymes", pd.Series([""] * len(t1))).map(
        lambda s: _E.get(str(s).strip(), ""))
    keep_map = {r["peptide"]: r for _, r in t1.iterrows()}
    t1_set = set(keep_map)

    # ---- 每肽聚合
    g = d.groupby("peptide").agg(
        n_spectra_isoPSV=("peptide", "size"),
        n_credible=("credible", "sum"),
        enzymes_isoPSV=("enzyme", lambda s: "+".join(sorted(set(s)))),
        novel_status=("novel_status", lambda s: ";".join(
            sorted({str(x) for x in s.dropna() if str(x)})[:2])),
    ).reset_index()
    g["n_credible"] = g["n_credible"].astype(int)
    g["in_list498"] = g["peptide"].map(lambda p: p in t1_set)
    g["table1_tier"] = g["peptide"].map(lambda p: keep_map.get(p, {}).get("tier", ""))
    g["table1_final_keep"] = g["peptide"].map(
        lambda p: keep_map.get(p, {}).get("final_keep", ""))
    g["table1_db_class"] = g["peptide"].map(
        lambda p: keep_map.get(p, {}).get("db_class", ""))
    g["isoPSV_verdict"] = "doubted"
    g.loc[g["n_credible"] > 0, "isoPSV_verdict"] = "credible"
    g = g.sort_values(["in_list498", "n_spectra_isoPSV"],
                      ascending=[False, False])

    # ---- 三级
    L1_pep, L1_sp = len(g), int(g["n_spectra_isoPSV"].sum())
    g2 = g[g["isoPSV_verdict"] == "credible"]
    L2_pep = len(g2)
    # 谱数用"判 credible 的谱"（行级），不是"通过肽的全部谱"，否则漏斗会虚高
    L2_sp = int(d["credible"].sum())
    L3 = g[g["in_list498"]]
    L3_pep = len(L3)
    L3_sp = int(L3["n_spectra_isoPSV"].sum())

    # 498 清单在 isoPSV 里的判读（按谱）
    l3_rows = d[d["peptide"].isin(t1_set)]
    verdict_sp = Counter("credible" if c else "doubted"
                         for c in l3_rows["credible"])
    doubt_sp = Counter()
    for x in l3_rows["novel_status"].dropna():
        for part in str(x).split(";"):
            if part:
                doubt_sp[part] += 1

    out = os.path.realpath(args.out)
    os.makedirs(out, exist_ok=True)
    cols = ["peptide", "enzymes_isoPSV", "n_spectra_isoPSV", "n_credible",
            "isoPSV_verdict", "novel_status", "in_list498", "table1_tier",
            "table1_final_keep", "table1_db_class"]
    g[cols].to_csv(os.path.join(out, "peptide_level_stats.tsv"), sep="\t",
                   index=False, encoding="utf-8-sig")
    try:
        with pd.ExcelWriter(os.path.join(out, "peptide_level_stats.xlsx")) as w:
            g[cols].to_excel(w, sheet_name="peptides", index=False)
            pd.DataFrame([
                {"口径": "1 所有新肽段（识别）", "肽段": L1_pep, "谱图": L1_sp},
                {"口径": "2 过滤新肽段（isoPSV credible）", "肽段": L2_pep, "谱图": L2_sp},
                {"口径": "3 保留新肽段（498 清单）", "肽段": L3_pep, "谱图": L3_sp},
                {"口径": "4 排除（所有 − 保留）", "肽段": L1_pep - L3_pep,
                 "谱图": L1_sp - L3_sp},
            ]).to_excel(w, sheet_name="summary", index=False)
            per.to_excel(w, sheet_name="by_enzyme", index=False)
    except Exception:
        pass

    # 按酶：保留（498）/ 排除（非 498）/ 合计，以及 credible 交叉
    # ORF 口径：保留用 Table1 的 `orf` 列（人工清单的指派）；排除/合计用原始 pFind 表
    # `Proteins` 的第一个 ORF（与 Table1 生成方式一致）。
    per = []
    for tag in ("Trypsin", "R29", "ALL"):
        s = d if tag == "ALL" else d[d["enzyme"] == tag]
        k = s[s["peptide"].isin(t1_set)]
        e = s[~s["peptide"].isin(t1_set)]
        c = s[s["credible"]]
        tt = t1 if tag == "ALL" else t1[t1["_enz"] == tag]
        orf_k = tt["orf"].nunique()
        orf_e = e["orf"].replace("", pd.NA).nunique()
        orf_all = len(set(tt["orf"].dropna()) | set(e["orf"].dropna()) - {""})
        per.append({
            "酶": tag,
            "保留_ORF": orf_k, "保留_谱": len(k), "保留_肽": k["peptide"].nunique(),
            "排除_ORF": orf_e, "排除_谱": len(e), "排除_肽": e["peptide"].nunique(),
            "合计_ORF": orf_all, "合计_谱": len(s), "合计_肽": s["peptide"].nunique(),
            "credible_谱": int(s["credible"].sum()), "credible_肽": c["peptide"].nunique(),
            "credible保留_谱": int(k["credible"].sum()),
            "credible保留_肽": k.loc[k["credible"], "peptide"].nunique(),
            "credible排除_谱": int(e["credible"].sum()),
            "credible排除_肽": e.loc[e["credible"], "peptide"].nunique(),
        })
    per = pd.DataFrame(per)
    per.to_csv(os.path.join(out, "stats_by_enzyme.tsv"), sep="\t",
               index=False, encoding="utf-8-sig")

    def prow(r):
        return ("| %s | %d / %d / %d | %d / %d / %d | %d / %d / %d |"
                % (r["酶"],
                   r["保留_ORF"], r["保留_肽"], r["保留_谱"],
                   r["排除_ORF"], r["排除_肽"], r["排除_谱"],
                   r["合计_ORF"], r["合计_肽"], r["合计_谱"]))

    # ---- 三分：已注释（命中注释库）/ 保留新肽段（498）/ 排除新肽段
    blob0 = ""
    if args.annotation and os.path.exists(args.annotation):
        blob0 = "|".join(iln(s) for s in _read_faa(args.annotation).values())
    d["in_db"] = d["peptide"].map(lambda p: bool(blob0) and p in blob0)
    d["keep"] = d["peptide"].isin(t1_set)
    # 保留优先（用户口径：498 条的谱 = 保留新肽段）；其余按命中与否分已注释/排除新肽段
    d["cat"] = "排除新肽段"
    d.loc[d["in_db"], "cat"] = "已注释"
    d.loc[d["keep"], "cat"] = "保留新肽段"
    five = d[d["keep"] & d["in_db"]]
    n_five = (five["peptide"].nunique(), len(five))

    cats = ("已注释", "排除新肽段", "保留新肽段")
    cat_per = []
    for tag in ("Trypsin", "R29", "ALL"):
        s = d if tag == "ALL" else d[d["enzyme"] == tag]
        row = {"酶": tag}
        for c in cats:
            x = s[s["cat"] == c]
            row[c + "_ORF"] = x["orf"].replace("", pd.NA).nunique()
            row[c + "_肽"] = x["peptide"].nunique()
            row[c + "_谱"] = len(x)
        cat_per.append(row)
    cat_per = pd.DataFrame(cat_per)
    cat_per.to_csv(os.path.join(out, "stats_by_category.tsv"), sep="\t",
                   index=False, encoding="utf-8-sig")

    lines = ["# 原始搜库结果（isoPSV 流程）· 肽段层面统计", "",
             "来源：`result_Trypsin(2)` / `result_R29` 的 `2.novel_peptide_competition_report.tsv`"
             "（isoPSV 新肽竞争报告，一行 = 一张谱，`a.Sequence` = ORF 库提出的候选肽）",
             "",
             "已剔除只在污染/自切库命中的 **%d 张谱**（`CON_*`）。" % n_con_rows,
             "",
             "**注释判定**：把每条肽回贴到 `%s`（%d 条蛋白），I/L 归一做子串匹配；"
             "命中 = 已注释肽段，未命中 = 新肽段。"
             % (os.path.basename(args.annotation or "注释库"),
                len(_read_faa(args.annotation)) if args.annotation else 0),
             "**新肽二分**：在 498 条人工清单内 = 保留新肽段（1,282 谱）；清单外 = 排除新肽段。",
             "",
             "所有表格式为 `ORF 数 / 肽段数 / 谱图数`；ORF 口径：保留按 Table1 的 `orf` 列，"
             "其余按原始 pFind 表 `Proteins` 的第一个 ORF（与 Table1 同法）。", "",
             "## 一、主表：已注释 / 排除新肽段 / 保留新肽段（按酶）", "",
             "| 酶 | 已注释 | 排除新肽段 | 保留新肽段 | 合计 |", "|---|---|---|---|---|"]
    for _, r in cat_per.iterrows():
        # 合计用 ORF 并集（各类 ORF 会重叠，直接相加会虚高）
        s = d if r["酶"] == "ALL" else d[d["enzyme"] == r["酶"]]
        tot = (s["orf"].replace("", pd.NA).nunique(),
               r["已注释_肽"] + r["排除新肽段_肽"] + r["保留新肽段_肽"],
               r["已注释_谱"] + r["排除新肽段_谱"] + r["保留新肽段_谱"])
        lines.append("| %s | %d / %d / %d | %d / %d / %d | %d / %d / %d | %d / %d / %d |"
                     % (r["酶"],
                        r["已注释_ORF"], r["已注释_肽"], r["已注释_谱"],
                        r["排除新肽段_ORF"], r["排除新肽段_肽"], r["排除新肽段_谱"],
                        r["保留新肽段_ORF"], r["保留新肽段_肽"], r["保留新肽段_谱"],
                        tot[0], tot[1], tot[2]))
    lines += ["",
              "> 保留新肽段 = 498 条人工清单（含其中 %d 条序列能命中注释库的，共 %d 谱）。"
              "若严格按「命中即已注释」优先，则保留变为 %d 肽 / %d 谱、已注释 +%d 肽 / +%d 谱。"
              % (n_five[0], n_five[1],
                 len(t1_set) - n_five[0],
                 int(cat_per.loc[cat_per["酶"] == "ALL", "保留新肽段_谱"].iloc[0]) - n_five[1],
                 n_five[0], n_five[1]),
              ""]
    lines += ["## 二、新肽总量与三分核对", "",
              "| 项 | 肽段 | 谱图 |", "|---|---|---|",
              "| 全部候选（ORF 命中，去污染） | %d | %d |" % (L1_pep, L1_sp),
              "| 已注释 | %d | %d |" % (int(cat_per.loc[cat_per["酶"] == "ALL",
                                                       "已注释_肽"].iloc[0]),
                                        int(cat_per.loc[cat_per["酶"] == "ALL",
                                                        "已注释_谱"].iloc[0])),
              "| 新肽合计（保留 + 排除） | %d | %d |"
              % (int(cat_per.loc[cat_per["酶"] == "ALL", "保留新肽段_肽"].iloc[0])
                 + int(cat_per.loc[cat_per["酶"] == "ALL", "排除新肽段_肽"].iloc[0]),
                 int(cat_per.loc[cat_per["酶"] == "ALL", "保留新肽段_谱"].iloc[0])
                 + int(cat_per.loc[cat_per["酶"] == "ALL", "排除新肽段_谱"].iloc[0])),
              "",
              "## 三、isoPSV credible × 三分（按酶）", "",
              "| 酶 | credible 合计 | 其中保留 | 其中排除新肽段 |",
              "|---|---|---|---|"]
    for _, r in per.iterrows():
        lines.append("| %s | %d / %d | %d / %d | %d / %d |"
                     % (r["酶"], r["credible_肽"], r["credible_谱"],
                        r["credible保留_肽"], r["credible保留_谱"],
                        r["credible排除_肽"], r["credible排除_谱"]))
    lines += ["", "## 四、三级（了解口径演变，可略）", "",
              "| 层级 | 肽段 | 谱图 |", "|---|---|---|",
              "| 1 所有新肽段（识别） | %d | %d |" % (L1_pep, L1_sp),
              "| 2 过滤新肽段（isoPSV credible） | %d | %d |" % (L2_pep, L2_sp),
              "| 3 保留新肽段（498 清单） | %d | %d |" % (L3_pep, L3_sp),
              "| 4 排除（所有 − 保留） | %d | %d |"
              % (L1_pep - L3_pep, L1_sp - L3_sp),
              ""]
    lines += ["", "## 四、按酶（旧口径的三列）", "",
              "| 酶 | 所有新肽段 | 过滤后（credible） | 498 清单 |",
              "|---|---|---|---|"]
    for _, r in per.iterrows():
        lines.append("| %s | %d 谱 / %d 肽 | %d 谱 / %d 肽 | %d 谱 / %d 肽 |"
                     % (r["酶"], r["合计_谱"], r["合计_肽"],
                        r["credible_谱"], r["credible_肽"],
                        r["保留_谱"], r["保留_肽"]))
    lines += ["", "## 五、498 清单在 isoPSV 里的判读（按谱）", "",
              "| 判读 | 谱数 |", "|---|---|",
              "| credible（通过） | %d |" % verdict_sp.get("credible", 0),
              "| doubted（存疑） | %d |" % verdict_sp.get("doubted", 0),
              "", "存疑原因（一谱可多条）：", "",
              "| 原因 | 谱数 |", "|---|---|"]
    for k, v in doubt_sp.most_common(8):
        lines.append("| %s | %d |" % (k, v))
    lines += ["",
              "## 六、说明", "",
              "1. **三分互斥**：已注释 / 排除新肽段 / 保留新肽段，三者相加 = 全部候选；"
              "「保留」优先于「已注释」——498 清单里若有肽能命中注释库（本轮 %d 条 / %d 谱），"
              "仍计入保留，表下已注明。" % (n_five[0], n_five[1]),
              "2. isoPSV 的 `credible` 是独立判定，见第三节交叉表——"
              "498 里只有 %d 条肽（%d 谱）被判 credible，"
              "清单**不是**从 isoPSV 通过集里筛出来的。"
              % (len(set(L3["peptide"]) & set(g2["peptide"])), verdict_sp.get("credible", 0)),
              "3. 每肽明细见 `peptide_level_stats.tsv`：`n_spectra_isoPSV` 是该肽在竞争报告里的谱数，"
              "`n_credible` 是其中判 credible 的谱数，`in_list498` 标记是否在 498 清单里。",
              ""]

    # ---- 序列层面旁证已在「一、已注释肽段」表中给出，此处不再重复
    with open(os.path.join(out, "stats_summary.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")

    print("\n".join(lines))
    print("→ %s" % out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
