# -*- coding: utf-8 -*-
"""QC 模块（设计书 v1.2：紧跟 E1 执行）。

输入：00_input/raw_info.tsv（方法参数+扫描统计）+ 04_search_annot 各任务
pFind-Filtered.spectra（Specificity / Miss.Clv.Sites 列做酶切质控）。
输出（04a_qc/）：qc_per_task.tsv（每 raw×CV 一行全指标）+ qc_summary.md
（门槛警告：ID 率<5%、IT 触顶>95% 且 AGC Fill<0.1、酶切特异性<90%）。
"""
import csv
import os
from collections import Counter

from .i18n import L


def parse_spectra(path):
    """pFind-Filtered.spectra → 统计 dict（PSM 数、特异性分布、漏切分布、肽段数）。"""
    if not os.path.exists(path):
        return None
    n, peps = 0, set()
    spec_c, mc_c = Counter(), Counter()
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            n += 1
            peps.add(row.get("Sequence", ""))
            spec_c[row.get("Specificity", "")] += 1
            mc_c[row.get("Miss.Clv.Sites", "")] += 1
    fully = spec_c.get("3", 0) + spec_c.get("5", 0)  # pFind：3/5=两端全切（TL）口径近似
    mc0 = mc_c.get("0", 0)
    return dict(n_psm=n, n_peptide=len(peps),
                specificity_dist=dict(spec_c), mc_dist=dict(mc_c),
                pct_full_cleavage=round(100 * fully / n, 2) if n else 0,
                pct_mc0=round(100 * mc0 / n, 2) if n else 0)


def build_qc(raw_info_tsv, annot_root, out_tsv, out_md):
    """汇总逐任务 QC 表与摘要。返回警告列表。"""
    rows = list(csv.DictReader(open(raw_info_tsv, newline="", encoding="utf-8-sig"),
                               delimiter="\t"))
    # raw_file → 方法/扫描统计
    info = {r["raw_file"]: r for r in rows}
    out_rows, warns = [], []
    if not os.path.isdir(annot_root):
        return warns
    for sp in sorted(os.listdir(annot_root)):
        sp_dir = os.path.join(annot_root, sp)
        if not os.path.isdir(sp_dir):
            continue
        for task in sorted(os.listdir(sp_dir)):
            tdir = os.path.join(sp_dir, task)
            spectra = os.path.join(tdir, "pFind-Filtered.spectra")
            st = parse_spectra(spectra)
            if st is None:
                continue
            short = task.rsplit("_", 1)[0] if "_" in task else task
            # 反查 raw（short 命名：物种_样品_酶_馏分_CV）——用 e1_report 若在
            e1r = os.path.join(annot_root, "e1_report.tsv")
            task2 = task
            # 指标：ID 率需要该任务对应 MGF 谱数（分母）——从 split_report 反查不可靠，
            # 改用 raw_info 该 raw 的 n_ms2 按 CV 比例不可行；此处先用 PSM 绝对值+肽数，
            # ID 率在 e1_report 有 n_psm；分母由调用方（cli）补。
            out_rows.append(dict(task=task2, species=sp, n_psm=st["n_psm"],
                                 n_peptide=st["n_peptide"],
                                 pct_full_cleavage=st["pct_full_cleavage"],
                                 pct_mc0=st["pct_mc0"],
                                 specificity_dist=str(st["specificity_dist"]),
                                 mc_dist=str(st["mc_dist"])))
            if st["pct_full_cleavage"] < 90:
                warns.append(f"{task2}: " + L("酶切特异性(两端全切近似)") +
                             f" {st['pct_full_cleavage']}% < 90%")
    fields = ["task", "species", "n_psm", "n_peptide", "pct_full_cleavage",
              "pct_mc0", "specificity_dist", "mc_dist"]
    with open(out_tsv, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(out_rows)
    # 方法级与扫描级警告（raw_info 汇总）
    for r in rows:
        sat = r.get("it_sat_pct", "")
        fill = r.get("med_fill_ms2", "")
        try:
            if sat and float(sat) > 95 and fill and float(fill) < 0.1:
                warns.append(f"{r['raw_file'][-40:]}: IT 触顶 {sat}% 且 MS2 AGC Fill {fill} <0.1")
        except ValueError:
            pass
    with open(out_md, "w", encoding="utf-8") as fh:
        fh.write(L(f"# QC 摘要（E1 后）") + "\n\n"
                 f"- {L('任务数（raw×CV）')}: {len(out_rows)}\n"
                 f"- {L('警告数')}: {len(warns)}\n\n" + L("## 警告") + "\n")
        for wn in warns:
            fh.write(f"- {L(wn)}\n")
        if not warns:
            fh.write(L("- 无") + "\n")
    return warns
