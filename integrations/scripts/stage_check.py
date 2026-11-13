#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""ORFscape 阶段自检：把"有产物但数字为 0 / 少一层"这类静默失败挡在路中间。

    python integrations/scripts/stage_check.py --workdir orfscape_projects/<proj>            # 全部阶段
    python integrations/scripts/stage_check.py --workdir orfscape_projects/<proj> --stage P6 # 单阶段

判定：FAIL 必须停；WARN 写进实验记录（多数 WARN 是静默丢失的早期信号）。
"""
import argparse
import csv
import os
import sys
from collections import Counter
from pathlib import Path

STAGES = ["P0", "P1", "D0", "C", "E1", "QC", "E2", "P6", "P6b", "P7", "P8", "P9"]


class Report:
    def __init__(self):
        self.rows = []

    def add(self, stage, level, msg):
        self.rows.append((stage, level, msg))

    def ok(self, stage, msg):
        self.add(stage, "PASS", msg)

    def warn(self, stage, msg):
        self.add(stage, "WARN", msg)

    def fail(self, stage, msg):
        self.add(stage, "FAIL", msg)

    def print(self, only=None):
        cur = None
        for stage, level, msg in self.rows:
            if only and stage != only:
                continue
            if stage != cur:
                print("\n== %s ==" % stage)
                cur = stage
            print("  [%s] %s" % (level, msg))
        bad = sum(1 for s, l, _ in self.rows if l == "FAIL" and (not only or s == only))
        warn = sum(1 for s, l, _ in self.rows if l == "WARN" and (not only or s == only))
        print("\n合计：%d FAIL / %d WARN" % (bad, warn))
        return bad


# ---------------------------------------------------------------- 读取工具
def faa_count(path):
    if not path.is_file():
        return 0
    n = 0
    with path.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                n += 1
    return n


def tsv_rows(path):
    if not path.is_file():
        return 0
    with path.open(encoding="utf-8-sig", errors="replace") as fh:
        return max(sum(1 for _ in fh) - 1, 0)


def read_tsv(path):
    if not path.is_file():
        return [], []
    with path.open(encoding="utf-8-sig", errors="replace", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return (list(rows[0].keys()) if rows else []), rows


def mgf_files(mgf_dir):
    if not mgf_dir.is_dir():
        return []
    return sorted(mgf_dir.glob("*.mgf"))


def task_dirs(base):
    if not base.is_dir():
        return []
    return sorted(d for d in base.iterdir() if d.is_dir())


def psm_count(base):
    """统计 base/<task>/pFind-Filtered.spectra 的数据行总数。"""
    total = 0
    for task in task_dirs(base):
        for f in task.rglob("pFind-Filtered.spectra"):
            with f.open(encoding="utf-8", errors="replace") as fh:
                total += max(sum(1 for _ in fh) - 1, 0)
    return total


def species_dirs(base):
    """物种子目录；跳过隐藏目录与下划线目录（如扫描器留下的 .mimosa）。"""
    if not base.is_dir():
        return []
    return sorted(d for d in base.iterdir()
                  if d.is_dir() and not d.name.startswith((".", "_")))


# ---------------------------------------------------------------- 各阶段
def check_p0(wd, rpt):
    info = wd / "00_input" / "raw_info.tsv"
    n = tsv_rows(info)
    if n > 0:
        rpt.ok("P0", "raw_info.tsv %d 行" % n)
    elif info.is_file():
        rpt.fail("P0", "raw_info.tsv 存在但 0 行")
    else:
        rpt.warn("P0", "还没有 scan 产物")
    return n


def check_p1(wd, rpt):
    for sp in species_dirs(wd / "01_genome"):
        miss = [f for f in ("genomic.fna", "genomic.gff", "protein.faa")
                if not (sp / f).is_file()]
        if miss:
            rpt.fail("P1", "%s 缺 %s" % (sp.name, ", ".join(miss)))
        else:
            rpt.ok("P1", "%s 三件套齐备" % sp.name)


def check_d0(wd, rpt, n_raw):
    mgfs = mgf_files(wd / "03_mgf")
    if not mgfs:
        rpt.fail("D0", "03_mgf 下没有 MGF")
        return 0
    rpt.ok("D0", "MGF %d 个" % len(mgfs))
    if n_raw and len(mgfs) != n_raw:
        rpt.fail("D0", "MGF 数 %d ≠ raw 数 %d" % (len(mgfs), n_raw))
    return len(mgfs)


def check_c(wd, rpt):
    # 只有带核心库（db_aug / L0）的目录才是"物种库"；只放辅助层（例如单独建的
    # L6 转录层）的子目录按辅助目录处理，不参与完整性检查（否则会误报 FAIL）。
    core, aux = [], []
    for sp in species_dirs(wd / "02_db"):
        names = {p.name for p in sp.iterdir() if p.is_file()}
        (core if ({"db_aug.faa", "L0_annot.faa"} & names) else aux).append(sp)
    if not core:
        rpt.fail("C", "02_db 下没有物种库（需要 db_aug.faa 或 L0_annot.faa）")
    for sp in aux:
        n = len([p for p in sp.glob("*.faa") if "decoy" not in p.name])
        rpt.ok("C", "%s 为辅助层目录（%d 个 FASTA），不做完整性检查" % (sp.name, n))
    for sp in core:
        d = sp
        n_l0 = faa_count(d / "L0_annot.faa")
        n_l2 = faa_count(d / "L2_6frame.faa")
        if n_l0 and n_l2:
            rpt.ok("C", "%s L0 %d / L2 %d" % (sp.name, n_l0, n_l2))
        else:
            rpt.fail("C", "%s L0=%d L2=%d（两层都必须非空）" % (sp.name, n_l0, n_l2))
        for layer, pat in (("L1", "L1_*.faa"), ("L3", "L3_altstart.faa"), ("L6", "L6_tx_orfs.faa")):
            n = sum(faa_count(p) for p in d.glob(pat) if "decoy" not in p.name)
            if layer == "L1" and n and not (d / "L1_coords.tsv").is_file():
                rpt.fail("C", "%s 有 L1 但缺 L1_coords.tsv（L1 无法定位）" % sp.name)
            elif n:
                rpt.ok("C", "%s %s %d 条" % (sp.name, layer, n))
        for lib in ("db_annot.faa", "db_aug.faa"):
            if not faa_count(d / lib):
                rpt.fail("C", "%s %s 为空或缺失" % (sp.name, lib))
        pools = [p for p in d.glob("db_6f*.faa") if "decoy" not in p.name]
        if not pools:
            rpt.fail("C", "%s 没有 db_6f 库" % sp.name)
        else:
            rpt.ok("C", "%s 六框库 %d 池" % (sp.name, len(pools)))


def check_e1(wd, rpt, n_mgf):
    base = wd / "04_search_annot"
    for sp in species_dirs(base):
        tasks = task_dirs(sp)
        psm = psm_count(sp)
        if not tasks:
            rpt.fail("E1", "%s 没有任务目录" % sp.name)
            continue
        rpt.ok("E1", "%s 任务 %d / PSM %d" % (sp.name, len(tasks), psm))
        if n_mgf and len(tasks) != n_mgf:
            rpt.fail("E1", "%s 任务数 %d ≠ MGF 数 %d" % (sp.name, len(tasks), n_mgf))
        if psm == 0:
            rpt.fail("E1", "%s PSM 为 0" % sp.name)


def check_qc(wd, rpt):
    qc = wd / "04a_qc"
    if not qc.is_dir() or not any(qc.iterdir()):
        rpt.warn("QC", "04a_qc 为空（qc / qcfigs / e1figs 未跑）")
    else:
        rpt.ok("QC", "04a_qc %d 个文件" % len(list(qc.iterdir())))


def check_e2(wd, rpt, n_mgf):
    base = wd / "04_search_disc"
    for sp in species_dirs(base):
        tasks = task_dirs(sp)
        psm = psm_count(sp)
        if not tasks:
            rpt.fail("E2", "%s 没有任务目录" % sp.name)
            continue
        rpt.ok("E2", "%s 任务 %d / PSM %d" % (sp.name, len(tasks), psm))
        if psm == 0:
            rpt.fail("E2", "%s PSM 为 0" % sp.name)
        # 库数从任务名尾部推断（形如 ...__6f1 / ...__aug），比数 faa 文件靠谱
        tags = set()
        for t in tasks:
            tags.add(t.name.split("__")[-1] if "__" in t.name else t.name)
        if n_mgf and len(tasks) % n_mgf != 0:
            rpt.warn("E2", "%s 任务数 %d 不能被 MGF 数 %d 整除" % (sp.name, len(tasks), n_mgf))
        elif n_mgf:
            rpt.ok("E2", "%s 每 raw 平均 %.1f 个库任务（%d 个库标签）"
                   % (sp.name, len(tasks) / n_mgf, len(tags)))


def check_p6(wd, rpt):
    for sp in species_dirs(wd / "05_novel"):
        p = sp / "novel_peptides.tsv"
        n = tsv_rows(p)
        if not n:
            rpt.fail("P6", "%s novel_peptides.tsv 为 0 行" % sp.name)
            continue
        _cols, rows = read_tsv(p)
        cnt = Counter()
        for r in rows:
            prot = str(r.get("protein") or "")
            for tag, pref in (("L1", "L1_"), ("L2", "L2_"), ("L3", "L3|"), ("L6", "L6|")):
                if prot.startswith(pref):
                    cnt[tag] += 1
                    break
        rpt.ok("P6", "%s 新肽 %d 行（层分布 %s）"
               % (sp.name, n, dict(cnt)))
        d = wd / "02_db" / sp.name
        for tag, pat in (("L1", "L1_*.faa"), ("L3", "L3_altstart.faa"), ("L6", "L6_tx_orfs.faa")):
            have = sum(faa_count(q) for q in d.glob(pat) if "decoy" not in q.name)
            if have and cnt.get(tag, 0) == 0:
                rpt.fail("P6", "%s %s 层有 %d 条 ORF 但 0 命中——白名单/改名 bug"
                         % (sp.name, tag, have))


def check_p6b(wd, rpt):
    for sp in species_dirs(wd / "05b_research"):
        p = sp / "novel_after_research.tsv"
        n = tsv_rows(p)
        p6 = tsv_rows(wd / "05_novel" / sp.name / "novel_peptides.tsv")
        if not n:
            rpt.warn("P6b", "%s 没有合并集（未跑 research 或全部失败）" % sp.name)
            continue
        if p6 and n < p6:
            rpt.fail("P6b", "%s 合并集 %d < P6 %d（回搜把肽弄丢了）" % (sp.name, n, p6))
        else:
            rpt.ok("P6b", "%s 合并集 %d 行（P6 %d）" % (sp.name, n, p6))


def check_p7(wd, rpt):
    p = wd / "06b_validate_research" / "novel_validated.tsv"
    if not p.is_file():
        p = wd / "06_validate" / "novel_validated.tsv"
    n = tsv_rows(p)
    if not n:
        rpt.fail("P7", "%s 不存在或 0 行" % p.parent.name)
        return
    _cols, rows = read_tsv(p)
    keep = sum(1 for r in rows if str(r.get("final_keep", "")).startswith("KEEP"))
    miss = sum(1 for r in rows if str(r.get("final_keep", "")).upper().startswith("MISS"))
    rpt.ok("P7", "%s：%d 行 / KEEP %d / MISSING %d" % (p.parent.name, n, keep, miss))
    if keep == 0:
        rpt.fail("P7", "KEEP 为 0——卡值或 isoPSV 环节有问题")
    if miss and "v_miss_reason" in _cols:
        no_reason = sum(1 for r in rows
                        if str(r.get("final_keep", "")).upper().startswith("MISS")
                        and not str(r.get("v_miss_reason") or "").strip())
        if no_reason:
            rpt.warn("P7", "%d 条 MISSING 没有原因" % no_reason)
    summary = p.parent / "validate_summary.md"
    if summary.is_file():
        for line in summary.read_text(encoding="utf-8", errors="replace").splitlines():
            if "卡值" in line:
                rpt.ok("P7", line.strip()[:110])
                break


def check_p8(wd, rpt):
    ng = wd / "07_novgene" / "NG_table.tsv"
    n = tsv_rows(ng)
    if not n:
        rpt.warn("P8", "还没有 07_novgene/NG_table.tsv")
    else:
        rpt.ok("P8", "NG_table %d 行（事件汇总）" % n)
    # 逐肽事件表在结果包里；这才是"未定位"会出现的地方
    per = wd / "P13_results" / "P13_3_events_per_peptide.tsv"
    if per.is_file():
        _cols, rows = read_tsv(per)
        bad = sum(1 for r in rows if any("未定位" in str(v) for v in r.values()))
        if bad:
            rpt.fail("P8", "逐肽事件表有 %d 行'未定位'（位点扫描漏了 05b_research）" % bad)
        else:
            rpt.ok("P8", "逐肽事件表 %d 行，无未定位" % len(rows))
    ann = wd / "08_annotation"
    if ann.is_dir() and any(ann.glob("*.tsv")):
        rpt.ok("P8", "08_annotation %d 个表" % len(list(ann.glob("*.tsv"))))
    else:
        rpt.warn("P8", "08_annotation 还没有表")


def check_p9(wd, rpt):
    rep = wd / "08_report"
    docs = list(rep.glob("*.docx")) if rep.is_dir() else []
    figs = list(rep.glob("*.png")) if rep.is_dir() else []
    if docs:
        rpt.ok("P9", "报告 %d 份、图 %d 张" % (len(docs), len(figs)))
    else:
        rpt.warn("P9", "08_report 还没有 docx（novreport 未跑）")


def main(argv=None):
    ap = argparse.ArgumentParser(description="ORFscape 阶段自检")
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--stage", default=None, choices=STAGES, help="只查某个阶段")
    args = ap.parse_args(argv)
    wd = Path(args.workdir).resolve()
    if not wd.is_dir():
        print("项目目录不存在：%s" % wd)
        return 2

    rpt = Report()
    n_raw = tsv_rows(wd / "00_input" / "raw_info.tsv")
    n_mgf = len(mgf_files(wd / "03_mgf"))
    only = args.stage

    def want(s):
        return only is None or only == s

    if want("P0"):
        check_p0(wd, rpt)
    if want("P1"):
        check_p1(wd, rpt)
    if want("D0"):
        check_d0(wd, rpt, n_raw)
    if want("C"):
        check_c(wd, rpt)
    if want("E1"):
        check_e1(wd, rpt, n_mgf)
    if want("QC"):
        check_qc(wd, rpt)
    if want("E2"):
        check_e2(wd, rpt, n_mgf)
    if want("P6"):
        check_p6(wd, rpt)
    if want("P6b"):
        check_p6b(wd, rpt)
    if want("P7"):
        check_p7(wd, rpt)
    if want("P8"):
        check_p8(wd, rpt)
    if want("P9"):
        check_p9(wd, rpt)

    print("ORFscape 阶段自检：%s" % wd)
    bad = rpt.print(only)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
