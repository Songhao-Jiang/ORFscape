# -*- coding: utf-8 -*-
"""把一次流程的用户侧结果打包到 `result_YYYYMMDD/`（YYYYMMDD = 任务开始日）。

只做**复制**，不改任何原始产物；目录内附 `README.md` 索引（数字、来源文件、
生成时间），拿到文件夹的人不需要进项目目录也能看懂。

    python -m orfscape pack --workdir <proj> [--out result_YYYYMMDD] [--date 20261114]

默认日期取 `run_state.json` 的 `created`（首次建项目的时间）；没有就看
`00_logs/` 里最早的日志时间；再不行用今天。
"""
import csv
import json
import shutil
import time
from pathlib import Path

# 用户侧最需要的几类产物（相对项目根）
TABLE_ITEMS = [
    ("06b_validate_research/novel_validated.tsv", "tables/novel_validated.tsv"),
    ("07_novgene/NG_table.tsv", "tables/NG_table.tsv"),
    ("07_novgene/NG_table.xlsx", "tables/NG_table.xlsx"),
    # 真核：事件表是 euk_events.tsv（七类、肽级归一），没有 NG_table
    ("07_novgene/euk_events.tsv", "tables/euk_events.tsv"),
    ("07_novgene/euk_events_summary.tsv", "tables/euk_events_summary.tsv"),
    ("08_annotation/AN_table.tsv", "tables/AN_table.tsv"),
    ("08_annotation/AN_table.xlsx", "tables/AN_table.xlsx"),
]


def _rows(path):
    with open(path, encoding="utf-8-sig", errors="replace", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def run_date(root):
    """任务开始日：run_state.created → 最早日志 → 今天。"""
    st = Path(root, "run_state.json")
    if st.is_file():
        try:
            d = json.loads(st.read_text(encoding="utf-8", errors="replace"))
            ts = (d.get("created") or "").strip()
            if len(ts) >= 10:
                return ts[:10].replace("-", "")
        except (OSError, ValueError):
            pass
    logs = sorted(Path(root, "00_logs").glob("*.log"), key=lambda p: p.stat().st_mtime)
    if logs:
        return time.strftime("%Y%m%d", time.localtime(logs[0].stat().st_mtime))
    return time.strftime("%Y%m%d")


def counts(root):
    """README 用的一页数字（能取到就取，取不到留空）。"""
    out = {}
    v = None
    for sub in ("06c_validate_full", "06b_validate_research", "06_validate"):
        cand = Path(root, sub, "novel_validated.tsv")
        if cand.is_file():
            v = cand
            break
    if v is not None:
        rows = _rows(v)
        keep = [r for r in rows if str(r.get("final_keep", "")).startswith("KEEP")]
        out["rows"] = len(rows)
        out["keep"] = len(keep)
        out["keep_unique"] = len({(r.get("peptide") or "").upper().replace("I", "L")
                                  for r in keep})
    ng = Path(root, "07_novgene", "NG_table.tsv")
    if ng.is_file():
        ev = {}
        for r in _rows(ng):
            name = (r.get("Event type") or "").strip()
            if not name:
                continue
            try:
                pep = int(float(str(r.get("Total pep") or 0).replace(",", "")))
                prot = int(float(str(r.get("Total prot") or 0).replace(",", "")))
            except ValueError:
                continue
            ev[name] = (pep, prot)
        out["events"] = ev
    if "events" not in out:
        euk = Path(root, "07_novgene", "euk_events.tsv")
        if euk.is_file():
            pep_c, orf_c = {}, {}
            for r in _rows(euk):
                name = (r.get("class") or "").strip()
                if not name:
                    continue
                pep_c[name] = pep_c.get(name, 0) + 1
                orf = (r.get("orf") or "").strip()
                orf_c.setdefault(name, set()).add(orf)
            out["events"] = {k: (v, len(orf_c.get(k, ()))) for k, v in pep_c.items()}
    an = Path(root, "08_annotation", "AN_table.tsv")
    if an.is_file():
        out["proteins"] = sum(1 for _ in an.open(encoding="utf-8", errors="replace")) - 1
    return out


def _en_summary(root, c, date):
    """English replacement for the Chinese validate_summary.md."""
    lines = ["# P7 validation summary (English)", "",
             "Run date: %s" % date,
             "Peptides (rows): %d; KEEP rows: %d; unique KEEP peptides: %d"
             % (c.get("rows", 0), c.get("keep", 0), c.get("keep_unique", 0)), "",
             "Retention rule: KEEP = final_keep; Tier1/Tier2 kept unconditionally, Tier3 must pass",
             "score / cosine / credible-novel gates. Cards are Tier1-anchored and not comparable",
             "across runs.", ""]
    for name in ("Intergenic new coding gene",
                 "N-term extension (same frame as annotated CDS)"):
        if c.get("events") and name in c["events"]:
            pep, prot = c["events"][name]
            lines.append("- %s: %d peptide records / %d ORFs" % (name, pep, prot))
    lines.append("")
    return chr(10).join(lines)


def write_readme(target, root, date, lang="en", log=print):
    c = counts(root)
    if lang == "zh":
        return _write_readme_zh(target, root, date, c, log)
    lines = ["# Results (result_%s)" % date, "",
             "**Project**: `%s`" % root,
             "**Packaged**: %s" % time.strftime("%Y-%m-%d %H:%M"),
             "**Run start**: %s" % date, "",
             "## Key numbers", "", "| Item | Value |", "|---|---|"]
    if "keep" in c:
        lines.append("| Candidate peptides (P7 input, rows) | %d |" % c.get("rows", 0))
        lines.append("| KEEP | %d rows / %d unique peptides |"
                     % (c["keep"], c.get("keep_unique", 0)))
    for name in ("Intergenic new coding gene",
                 "N-term extension (same frame as annotated CDS)"):
        if c.get("events") and name in c["events"]:
            pep, prot = c["events"][name]
            lines.append("| %s | %d peptides / %d ORFs |" % (name, pep, prot))
    if "proteins" in c:
        lines.append("| Annotated novel proteins | %d |" % c["proteins"])
    lines += ["", "## Contents", "",
              "| Path | Content |", "|---|---|",
              "| `NTM_report_EN.md` / `.docx` | Full pipeline report (figures, legends, analysis) |",
              "| `NovelGene_Report.docx` | Auto-generated summary of the novel-gene analysis |",
              "| `tables/` | Validated peptides, event table, annotation table, validation summary |",
              "| `figures/` | Report figures (PNG) |",
              "| `pdf/` | Peptide/IGV booklets if this project produced them |",
              "| `README.md` | This file |", "",
              "## Conventions", "",
              "- KEEP = P7 `final_keep`; unique peptides (I/L-normalised) are reported externally,",
              "  row counts are used for ORF-level events.",
              "- Events are exclusive per ORF (N-terminal correction > frame > antisense > intergenic).",
              "- Card values are Tier1-anchored and **not comparable across runs**.", ""]
    Path(target, "README.md").write_text(chr(10).join(lines), encoding="utf-8")
    log("  README → %s" % Path(target, "README.md"))


def _write_readme_zh(target, root, date, c, log=print):
    lines = ["# 结果说明（result_%s）" % date, "",
             "**项目**：`%s`" % root,
             "**打包时间**：%s" % time.strftime("%Y-%m-%d %H:%M"),
             "**任务开始**：%s" % date, "",
             "## 一、关键数字", "", "| 项 | 值 |", "|---|---|"]
    if "keep" in c:
        lines.append("| 候选肽（P7 输入） | %d |" % c.get("rows", 0))
        lines.append("| 保留（KEEP） | %d 行 / %d 唯一肽 |"
                     % (c["keep"], c.get("keep_unique", 0)))
    for name in ("Intergenic new coding gene",
                 "N-term extension (same frame as annotated CDS)"):
        if c.get("events") and name in c["events"]:
            pep, prot = c["events"][name]
            lines.append("| %s | %d 肽 / %d ORF |" % (name, pep, prot))
    if "proteins" in c:
        lines.append("| 新编码蛋白注释 | %d |" % c["proteins"])
    lines += ["", "## 二、目录内容", "",
              "| 目录 | 内容 |", "|---|---|",
              "| `report/` | Word 汇总报告 |",
              "| `tables/` | 逐肽结论、事件表、注释表、收口汇总 |",
              "| `figures/` | 汇报用统计图（PNG） |",
              "| `pdf/` | 谱册 / 事件册（若本项目生成过） |",
              "| `README.md` | 本说明 |", "",
              "## 三、口径提醒", "",
              "- 保留 = P7 `final_keep=KEEP`；对外计**唯一肽**（I/L 归一）。",
              "- 事件按 ORF 互斥归类（N 端校正 > 不同框 > 反义 > 基因间区）。",
              "- 卡值按本轮 Tier1 标定，**跨轮不可比**；换库/换版必须整链重跑。", ""]
    p = Path(target, "README.md")
    p.write_text("\n".join(lines), encoding="utf-8")
    log("  README → %s" % p)


def pack(root, out=None, date=None, lang="en", log=print):
    root = Path(root).resolve()
    if not root.is_dir():
        raise SystemExit("项目目录不存在：%s" % root)
    date = date or run_date(root)
    target = Path(out) if out else Path(root, "result_%s" % date)
    target = target.resolve()
    if target == root or root not in target.parents:
        raise SystemExit("输出目录必须位于项目内：%s" % target)
    for sub in ("report", "tables", "figures", "pdf"):
        Path(target, sub).mkdir(parents=True, exist_ok=True)

    copied = []
    for rel, dst in TABLE_ITEMS:
        src = Path(root, rel)
        if src.is_file():
            shutil.copy2(src, Path(target, dst))
            copied.append(dst)
    if lang == "zh":
        ngr = Path(root, "08_report_zh", "NovelGene_Report_zh.docx")
        if ngr.is_file():
            shutil.copy2(ngr, Path(target, "report", "NovelGene_Report_zh.docx"))
            copied.append("report/NovelGene_Report_zh.docx")
    else:
        ngr = Path(root, "08_report", "NovelGene_Report.docx")
        if ngr.is_file():
            shutil.copy2(ngr, Path(target, "report", "NovelGene_Report.docx"))
            copied.append("report/NovelGene_Report.docx")
    if lang == "zh":
        vs = Path(root, "06b_validate_research", "validate_summary.md")
        if vs.is_file():
            shutil.copy2(vs, Path(target, "tables", "validate_summary.md"))
            copied.append("tables/validate_summary.md")
    else:
        Path(target, "tables", "validate_summary_EN.md").write_text(
            _en_summary(root, counts(root), date), encoding="utf-8")
        copied.append("tables/validate_summary_EN.md")

    # 图：只收本项目的 docs/figs_<前缀>（不混别的项目，也不收 08_report 的历史图）
    prefix = root.name.split("_")[0]
    fig_dir = Path(root).parent.parent / "ORFscape" / "docs" / ("figs_" + prefix.lower()
                                                       + ("_zh" if lang == "zh" else ""))
    if fig_dir.is_dir():
        for png in sorted(fig_dir.glob("*.png")):
            shutil.copy2(png, Path(target, "figures", png.name))
            copied.append("figures/" + png.name)
    # 事件/注释系列图（novgene / novannot 产出，报告里会引用）
    ng_dir, an_dir = ("07_novgene_zh", "08_annotation_zh") if lang == "zh"         else ("07_novgene", "08_annotation")
    for rel in (ng_dir + "/NG1_events_bar.png", ng_dir + "/NG2_orf_distributions.png",
                an_dir + "/AN1_properties.png", an_dir + "/AN2_homology.png",
                an_dir + "/AN3_features.png"):
        src = Path(root, rel)
        if src.is_file():
            shutil.copy2(src, Path(target, "figures", src.name))
            copied.append("figures/" + src.name)

    # PDF（TP 这类有册子的项目）
    for pdf in sorted(Path(root, "P13_results").glob("*.pdf")):
        shutil.copy2(pdf, Path(target, "pdf", pdf.name))
        copied.append("pdf/" + pdf.name)

    # 报告：按语言选（EN 默认），以固定英文名放到包根，图用 figures/ 相对路径
    prefix_up = prefix.upper()
    tag = "ZH" if lang == "zh" else "EN"
    # 只拷 md：docx 要在包内用相对 figures/ 路径重建（图片才能真正嵌进去）
    pats = ("%s_*全流程报告_%s_*.md" % (prefix_up, tag),)
    got = False
    for pat in pats:
        cands = sorted(Path(root).parent.parent.glob(pat), key=lambda q: q.stat().st_mtime)
        for src in cands:
            dst = Path(target, "%s_report_%s%s" % (prefix_up, tag, src.suffix))
            shutil.copy2(src, dst)
            copied.append(dst.name)
            got = True
    if not got:
        for pat in ("%s_*全流程报告*.md" % prefix_up,):
            cands = sorted(Path(root).parent.parent.glob(pat),
                           key=lambda q: q.stat().st_mtime)
            if cands:
                src = cands[-1]
                dst = Path(target, "%s_report%s" % (prefix_up, src.suffix))
                shutil.copy2(src, dst)
                copied.append(dst.name)

    write_readme(target, root, date, lang=lang, log=log)
    log("  %d 个文件 → %s" % (len(copied), target))
    return str(target)
