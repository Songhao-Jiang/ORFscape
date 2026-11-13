# -*- coding: utf-8 -*-
"""阶段产物自检（GUI/编排器内置版）。

定位：**便宜、可重复、抓静默失败**。检查的是"文件在不在、里面有东西没有、
数量级对不对"，不做生物学判读；深度判读仍由 skill 的 ``stage_check.py`` 负责。

    python -m orfscape.checks --workdir <proj> --stage C
    python -m orfscape.checks --workdir <proj>            # 全部阶段

返回：``{"status": ok|warn|fail|skip, "fails": n, "warns": n, "lines": [...]}``
"""
import argparse
import glob
import json
import os
import sys

PARENT = os.pardir


def _root(workdir):
    raw = str(workdir).replace("\\", "/")
    if any(part == PARENT for part in raw.split("/")):
        raise ValueError("project path contains a parent component")
    root = os.path.realpath(os.path.abspath(str(workdir)))
    if not os.path.isdir(root):
        raise ValueError("project directory not found")
    return root


def _n_lines(path):
    try:
        with open(path, "rb") as fh:
            return sum(1 for _ in fh)
    except OSError:
        return 0


def _rows(path):
    """表格数据行数（扣表头）。"""
    return max(0, _n_lines(path) - 1)


def _faa_count(path):
    """FASTA 记录数（> 开头行）。"""
    n = 0
    try:
        with open(path, "rb") as fh:
            for line in fh:
                if line.startswith(b">"):
                    n += 1
    except OSError:
        return 0
    return n


def _size(path):
    try:
        return os.path.getsize(path)
    except OSError:
        return 0


def _glob(root, pattern):
    base = os.path.realpath(root)
    hits = []
    for p in glob.glob(os.path.join(base, pattern)):
        rp = os.path.realpath(p)
        if rp == base or rp.startswith(base + os.sep):
            hits.append(rp)
    return sorted(hits)


def _psm_total(files):
    from .core import pfsearch
    return sum(max(0, pfsearch.count_result_rows(f)) for f in files)


class Report(object):
    def __init__(self, stage):
        self.stage = stage
        self.lines = []
        self.fails = 0
        self.warns = 0

    def ok(self, text):
        self.lines.append("  ok   " + text)

    def warn(self, text):
        self.warns += 1
        self.lines.append("  WARN " + text)

    def fail(self, text):
        self.fails += 1
        self.lines.append("  FAIL " + text)

    def result(self):
        if self.fails:
            status = "fail"
        elif self.warns:
            status = "warn"
        else:
            status = "ok"
        summary = "%s: %d fail, %d warn" % (self.stage, self.fails, self.warns)
        return {"stage": self.stage, "status": status, "fails": self.fails,
                "warns": self.warns, "lines": self.lines, "summary": summary}


def _require(report, path, label, min_bytes=1):
    if not os.path.isfile(path):
        report.fail("%s missing" % label)
        return False
    if _size(path) < min_bytes:
        report.fail("%s is empty" % label)
        return False
    report.ok("%s" % label)
    return True


def check_stage(workdir, stage):
    """按阶段做断言；未知阶段返回 skip。"""
    try:
        root = _root(workdir)
    except ValueError as exc:
        return {"stage": stage, "status": "fail", "fails": 1, "warns": 0,
                "lines": ["  FAIL " + str(exc)], "summary": str(exc)}
    rp = Report(stage)

    if stage == "P0":
        _require(rp, os.path.join(root, "00_input", "raw_info.tsv"), "raw_info.tsv", 100)
        sp = os.path.join(root, "00_input", "species.tsv")
        _require(rp, sp, "species.tsv", 20)
        n_sp = max(0, _rows(sp))
        if n_sp == 0:
            rp.fail("no species resolved")
        grp = os.path.join(root, "00_input", "groups.csv")
        if os.path.isfile(grp):
            n = max(0, _rows(grp))
            if n == 0:
                rp.fail("groups.csv has no rows")
            else:
                rp.ok("%d grouped raw file(s)" % n)
        else:
            rp.warn("groups.csv missing")

    elif stage == "P1":
        dirs = [d for d in _glob(root, "01_genome/*") if os.path.isdir(d)]
        if not dirs:
            rp.fail("01_genome/<species> missing")
        for d in dirs:
            code = os.path.basename(d)
            _require(rp, os.path.join(d, "genomic.fna"), "%s genomic.fna" % code, 1000)
            _require(rp, os.path.join(d, "protein.faa"), "%s protein.faa" % code, 1000)
            if not os.path.isfile(os.path.join(d, "genomic.gff")):
                rp.warn("%s genomic.gff missing (alt-frame events will be weaker)" % code)

    elif stage == "D0":
        mgfs = _glob(root, "03_mgf/*.mgf")
        if not mgfs:
            rp.fail("no MGF in 03_mgf/")
        else:
            total = sum(_size(f) for f in mgfs) / (1024.0 ** 3)
            rp.ok("%d MGF, %.2f GB" % (len(mgfs), total))

    elif stage == "C":
        dbdirs = [d for d in _glob(root, "02_db/*") if os.path.isdir(d)]
        if not dbdirs:
            rp.fail("02_db/<species> missing")
        # 只有带核心库（db_aug/db_annot/L0）的目录才是"物种库"，
        # 其余（如只放新 L6 层的 TP_l6b）按辅助层目录处理，不判失败。
        core, aux = [], []
        for d in dbdirs:
            names = set(os.listdir(d))
            (core if ({"db_aug.faa", "L0_annot.faa"} & names) else aux).append(d)
        if dbdirs and not core:
            rp.fail("no species database under 02_db/ (needs db_aug.faa or L0_annot.faa)")
        for d in core:
            code = os.path.basename(d)
            for name, min_b in (("L0_annot.faa", 100), ("L2_6frame.faa", 100),
                                ("db_annot.faa", 100), ("db_aug.faa", 100),
                                ("db_6f.faa", 100)):
                _require(rp, os.path.join(d, name), "%s/%s" % (code, name), min_b)
            for name in ("L1_prodigal.faa", "L3_altstart.faa", "L6_tx_orfs.faa"):
                p = os.path.join(d, name)
                if not os.path.isfile(p):
                    rp.warn("%s/%s missing (layer not built)" % (code, name))
                elif _size(p) < 60:
                    rp.fail("%s/%s is empty — the layer was built but has no sequence"
                            % (code, name))
        for d in aux:
            code = os.path.basename(d)
            n_faa = len([f for f in os.listdir(d) if f.endswith(".faa")])
            if n_faa:
                rp.warn("%s/ is an auxiliary layer dir (%d FASTA, no core database)"
                        % (code, n_faa))

    elif stage in ("E1", "E2", "QC"):
        sub = {"E1": "04_search_annot", "E2": "04_search_disc", "QC": "04a_qc"}[stage]
        base = os.path.join(root, sub)
        if stage == "QC":
            files = _glob(root, "04a_qc/*")
            if not files:
                rp.warn("04a_qc/ is empty (qcfigs/e1figs not run)")
            else:
                rp.ok("%d QC output(s)" % len(files))
        else:
            specs = _glob(root, sub + "/*/*/pFind-Filtered.spectra")
            if not specs:
                rp.fail("no pFind-Filtered.spectra under %s/" % sub)
            else:
                n_psm = _psm_total(specs)
                if n_psm <= 0:
                    rp.fail("%d result file(s) but 0 PSM" % len(specs))
                else:
                    rp.ok("%d task(s), %d PSM" % (len(specs), n_psm))
                n_mgf = len(_glob(root, "03_mgf/*.mgf"))
                if n_mgf and len(specs) < n_mgf:
                    rp.warn("%d task(s) for %d MGF (partial run?)"
                            % (len(specs), n_mgf))
            if not os.path.isdir(base):
                rp.fail("%s/ missing" % sub)

    elif stage == "EUK":
        # 真核层：eukbuild 产物（L0_annot / L0_junction / L2_txorf + manifest）
        found = False
        for d in [x for x in _glob(root, "02_db/*") if os.path.isdir(x)]:
            mf = os.path.join(d, "euk_manifest.json")
            if not os.path.isfile(mf):
                continue
            found = True
            code = os.path.basename(d)
            try:
                with open(mf, encoding="utf-8") as fh:
                    man = json.load(fh)
            except (OSError, ValueError) as exc:
                rp.fail("%s euk_manifest.json 读不了：%s" % (code, exc))
                continue
            counts = man.get("counts") or {}
            n_prot = counts.get("proteins") or 0
            n_junc = counts.get("junction_peptides") or 0
            real_prot = _faa_count(os.path.join(d, "L0_annot.faa"))
            if n_prot <= 0 or real_prot <= 0:
                rp.fail("%s 注释蛋白为空（manifest=%d, 文件=%d）"
                        % (code, n_prot, real_prot))
            elif real_prot != n_prot:
                rp.warn("%s 蛋白数与 manifest 不一致（%d vs %d）"
                        % (code, real_prot, n_prot))
            else:
                rp.ok("%s 注释蛋白 %d 条" % (code, n_prot))
            if not os.path.isfile(os.path.join(d, "L0_junction.faa")):
                rp.fail("%s 缺 L0_junction.faa" % code)
            elif n_junc <= 0:
                rp.warn("%s 没有造出剪接连接肽（单外显子物种可忽略）" % code)
            else:
                rp.ok("%s 剪接连接肽 %d 条（唯一 %d）"
                      % (code, n_junc, counts.get("junction_unique") or 0))
            pre = (man.get("preset") or {}).get("check") or {}
            if pre.get("expected"):
                if pre.get("ok"):
                    rp.ok("%s 蛋白数与文档一致（%d vs %d，%+.2f%%）"
                          % (code, pre["observed"], pre["expected"],
                             pre["delta_pct"]))
                else:
                    rp.warn("%s 蛋白数与文档偏差 %+.2f%%（文档 %d）"
                            % (code, pre["delta_pct"], pre["expected"]))
            rc = man.get("reference_check") or {}
            if rc:
                rate = rc.get("exact_rate") or 0
                match = rc.get("match_rate") or 0
                if rate >= 0.99:
                    rp.ok("%s 与参考蛋白完全一致 %.2f%%" % (code, rate * 100))
                elif match >= 0.99:
                    # 人类这类：exact 低是 N 端起始 codon（RefSeq 强制 M）
                    # 与基因组-参考替换造成的，match_rate 已把它们算进来
                    rp.ok("%s 与参考蛋白一致 %.2f%%（完全一致 %.2f%%；"
                          "差异集中在 N 端起始与参考替换，非剪接错误）"
                          % (code, match * 100, rate * 100))
                elif rate >= 0.95:
                    rp.warn("%s 与参考蛋白一致 %.2f%%，含替换重建 %.2f%%"
                            "（检查密码子表/起始）"
                            % (code, rate * 100, match * 100))
                else:
                    rp.fail("%s 与参考蛋白完全一致只有 %.2f%%（剪接可能有错）"
                            % (code, rate * 100))
        if not found:
            rp.fail("没有真核建库产物（02_db/*/euk_manifest.json）")

    elif stage == "P6":
        tabs = _glob(root, "05_novel/*/novel_peptides.tsv")
        if not tabs:
            rp.fail("05_novel/*/novel_peptides.tsv missing")
        for t in tabs:
            n = _rows(t)
            if n <= 0:
                rp.fail("%s has 0 peptides" % os.path.basename(os.path.dirname(t)))
            else:
                rp.ok("%s: %d peptide(s)" % (os.path.basename(os.path.dirname(t)), n))

    elif stage == "P6b":
        tabs = _glob(root, "05b_research/*/novel_after_research.tsv")
        if not tabs:
            rp.fail("05b_research/*/novel_after_research.tsv missing")
        for t in tabs:
            n = _rows(t)
            (rp.ok if n > 0 else rp.fail)("%s: %d peptide(s)"
                                          % (os.path.basename(os.path.dirname(t)), n))

    elif stage == "P7":
        t = os.path.join(root, "06b_validate_research", "novel_validated.tsv")
        if not os.path.isfile(t):
            rp.fail("06b_validate_research/novel_validated.tsv missing")
        else:
            n = _rows(t)
            (rp.ok if n > 0 else rp.fail)("validated peptides: %d" % n)

    elif stage == "P8":
        for rel in ("07_novgene/NG_table.tsv",):
            p = os.path.join(root, *rel.split("/"))
            if not os.path.isfile(p):
                rp.fail(rel + " missing")
            else:
                n = _rows(p)
                (rp.ok if n > 0 else rp.fail)("%s: %d event(s)" % (rel, n))
        ann = _glob(root, "08_annotation/*.tsv")
        if not ann:
            rp.warn("08_annotation/*.tsv missing")
        else:
            rp.ok("%d annotation table(s)" % len(ann))

    elif stage == "P9":
        docs = _glob(root, "08_report/*.docx") + _glob(root, "08_report/*.pdf")
        if not docs:
            rp.warn("08_report/ has no docx/pdf yet")
        else:
            rp.ok("%d report file(s)" % len(docs))

    else:
        return {"stage": stage, "status": "skip", "fails": 0, "warns": 0,
                "lines": [], "summary": "%s: no checks defined" % stage}

    return rp.result()


def check_all(workdir, stages=None):
    stages = stages or ["P0", "P1", "D0", "C", "E1", "QC", "E2", "P6", "P6b",
                        "P7", "P8", "P9"]
    return [check_stage(workdir, s) for s in stages]


def main(argv=None):
    ap = argparse.ArgumentParser(description="ORFscape stage checks")
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--stage", default=None, help="single stage; default: all")
    ap.add_argument("--json", action="store_true", help="print JSON")
    args = ap.parse_args(argv)
    if args.stage:
        res = [check_stage(args.workdir, args.stage)]
    else:
        res = check_all(args.workdir)
    if args.json:
        print(json.dumps(res, ensure_ascii=False, indent=1))
    else:
        for r in res:
            print("[%s] %s" % (r["status"].upper(), r["summary"]))
            for line in r["lines"]:
                print(line)
    return 0 if all(r["status"] != "fail" for r in res) else 1


if __name__ == "__main__":
    sys.exit(main())
