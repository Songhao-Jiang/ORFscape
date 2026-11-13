# -*- coding: utf-8 -*-
"""SRA 转录组自动挑组（第一原则：全自动）。

策略：按 BioProject 分组（不同 study≈不同条件）→ 每 study 取 ≥3 run（配对/适中优先）
     → 多 study 合计 ≤5 组 × ≤4 run = ≤20 run。
"""
import csv
import os

from . import sraselect


def auto_select(taxid, out_tsv, max_studies=5, min_reps=3, max_reps=4):
    rows, by_bp = sraselect.fetch_runs(taxid)
    # 排序：study 内 run 数多的优先
    ranked = sorted(by_bp.items(), key=lambda kv: -len(kv[1]))
    selected = []
    for bp, runs in ranked:
        if len(runs) < min_reps:
            continue
        if len(selected) >= max_studies:
            break
        # 排序：配对优先 → 适中数据量（5-30 Gb）
        def sort_key(r):
            paired = 0 if r["layout"] == "PAIRED" else 1
            mb = r["mb"] or 0
            size_penalty = abs(mb - 10000)
            return (paired, size_penalty)
        picks = sorted(runs, key=sort_key)[:max_reps]
        selected.append((bp, picks))
    # 落盘
    out_rows = []
    for bp, picks in selected:
        for r in picks:
            out_rows.append([r["run"], bp, r["biosample"], r["layout"],
                             r["spots"], r["mb"], r["title"], r["instrument"]])
    with open(out_tsv, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["run", "bioproject", "biosample", "layout", "spots", "Mb",
                    "title", "instrument"])
        w.writerows(out_rows)
    return out_rows, selected


def summary(selected):
    lines = []
    for i, (bp, picks) in enumerate(selected, 1):
        t = picks[0]["title"][:50] or "n/a"
        lines.append(f"  组{i} [{bp}] {len(picks)} runs | {t}")
    return "\n".join(lines)
