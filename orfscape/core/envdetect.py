# -*- coding: utf-8 -*-
"""环境检测器（设计书 A 阶段）：外部工具/引擎的自动发现与版本探测。

设计要点：
- 探测两类来源：固定根目录 glob（本机 Bioinformatics 布局）+ PATH（shutil.which）；
- 只读探测，不安装（安装需用户在 GUI/CLI 显式允许，A 阶段弹窗逻辑在 GUI 轮实现）；
- 产出 env_report.tsv（人读）+ run_state 登记（机器读），缺失项给出"将自动安装"分级。

用法（cli: doctor 子命令）：
  rows = detect_all()                 # 每项 dict: name/category/path/version/status/note
  write_report(rows, out_tsv)
"""
import os
import re
import shutil
import subprocess

from .. import localpaths


def _roots():
    """探测根：ORFSCAPE_TOOL_ROOTS（os.pathsep 分隔）在前，其后常见安装位置。"""
    env = os.environ.get("ORFSCAPE_TOOL_ROOTS") or ""
    roots = [p for p in env.split(os.pathsep) if p]
    roots += [
        r"D:\Bioinformatics",
        r"C:\Program Files\ProteoWizard",
        r"C:\Program Files\NCBI",
        r"C:\software",
        r"D:\isoPSV_bin",
    ]
    iso = localpaths.get("isopsv")
    if iso:
        roots.append(iso)
    return roots


# 固定探测根（工具装在别处：设 ORFSCAPE_TOOL_ROOTS 或 paths.local.json）
_ROOTS = _roots()


def _glob_first(pattern):
    """在固定根下按 glob 找第一个匹配（不递归全盘）。"""
    import glob
    for root in _ROOTS:
        hits = sorted(glob.glob(os.path.join(root, pattern)))
        if hits:
            return hits[-1]  # 取字典序最大（通常=版本号最新）
    return ""


# ---------------------------------------------------------------- 探测项注册表
# kind: exe（可执行路径探测）| dir（目录存在即可）
# need: required（流程必需）/ recommended（仲裁/建库增强）/ optional（可选槽位）
PROBES = [
    dict(name="pFind", category="发现搜库", need="required", patterns=[
        r"pFind3.1\bin\Searcher.exe", r"pFind*\bin\Searcher.exe"]),
    dict(name="pParse", category="谱图转换", need="required", patterns=[
        r"pFind3.1\bin\pParse.exe", r"pFind*\bin\pParse.exe"]),
    dict(name="msconvert(ProteoWizard)", category="电压拆分/转换", need="required", patterns=[
        r"ProteoWizard\msconvert.exe", r"ProteoWizard*\msconvert.exe",
        r"ProteoWizard*\ProteoWizard*\bin\msconvert.exe"], which="msconvert"),
    dict(name="MaxQuant", category="免费引擎(仲裁槽)", need="optional", patterns=[
        r"MaxQuant*\MaxQuant*\bin\MaxQuantCmd.exe", r"MaxQuant*\bin\MaxQuantCmd.exe"], multi=True),
    dict(name="Comet", category="免费引擎(仲裁槽)", need="recommended", patterns=[
        r"Comet*\comet*.exe"], which="comet"),
    dict(name="MSFragger", category="免费引擎(仲裁槽)", need="recommended", patterns=[
        r"MSfragger\MSFragger-*\MSFragger*.jar", r"MSfragger*\MSFragger*.jar"]),
    dict(name="Sage", category="免费引擎(仲裁槽)", need="recommended", patterns=[
        r"Sage\sage-*\sage.exe", r"Sage*\sage.exe", r"Sage*\sage*.exe"], which="sage"),
    dict(name="MetaMorpheus", category="免费引擎(仲裁槽)", need="optional", patterns=[
        r"MetaMorpheus*\MetaMorpheus_CommandLine\CMD.exe", r"MetaMorpheus*\MetaMorpheus*\MetaMorpheus.exe"]),
    dict(name="AlphaPeptDeep", category="谱图预测器", need="recommended", patterns=[
        r"AlphaPeptDeep*\peptdeep*.exe"]),
    dict(name="BLAST+", category="功能注释", need="recommended", patterns=[
        r"BLAST\blastp.exe", r"blast*\blastp.exe"], which="blastp"),
    dict(name="Prodigal", category="ab initio 预测(L1)", need="recommended", patterns=[
        r"Prodigal*\prodigal*.exe"], which="prodigal"),
    dict(name="Augustus", category="ab initio 预测(L1)", need="optional", patterns=[
        r"augustus*\bin\augustus.exe", r"augustus*\augustus.exe"], which="augustus"),
    dict(name="GeneMarkS-2", category="ab initio 预测(L1)", need="optional", patterns=[
        r"gmst*\gmst.exe", r"GeneMarkS*\gmst*.exe"]),
    dict(name="rnaSPAdes", category="转录本层(L6)", need="optional", patterns=[
        r"SPAdes*\bin\spades.exe"], which="spades"),
    dict(name="isoPSV", category="竞争仲裁后端", need="required", patterns=[
        r"isoPSV.exe", r"isoPSV\isoPSV.exe", r"isoPSV\start_isoPSV.bat", r"isoPSV*\start_isoPSV.bat"]),
    dict(name="ThermoRawFileParser", category="谱图转换回退", need="optional", patterns=[
        r"ThermoRawFileParser*\ThermoRawFileParser.exe"], which="ThermoRawFileParser"),
    dict(name="rawmethod(vendored)", category="QC方法提取", need="required", patterns=[], vendored=True),
    dict(name="rawdump(vendored)", category="QC扫描表", need="required", patterns=[], vendored=True),
]

# 降级链说明（检测缺失时提示影响）
DEGRADE = {
    "Prodigal": "L1 层降级：无 Prodigal 则跳过 ab initio 层（六框层部分覆盖）",
    "Augustus": "L1 层降级：仅 Prodigal（如 Prodigal 也在缺则 L1 整层跳过）",
    "GeneMarkS-2": "L1 层降级：仅 Prodigal+Augustus",
    "rnaSPAdes": "L6 转录本层跳过（报告声明）",
    "MaxQuant": "仲裁引擎槽减一（Comet/Sage/MSFragger 可用即不受阻）",
    "Comet": "仲裁引擎槽减一",
    "MSFragger": "仲裁引擎槽减一",
    "Sage": "仲裁引擎槽减一",
    "AlphaPeptDeep": "预测器回退 pDeep2（如配）或无预测模式",
    "BLAST+": "功能注释降级（跳过 BLASTp 层）",
}


def _vendored_path(name):
    base = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # orfscape/
    if "rawmethod" in name:
        p = os.path.join(base, "tools", "rawmethod", "rawmethod.exe")
    else:
        p = os.path.join(base, "tools", "rawdump", "rawdump.exe")
    return p if os.path.isfile(p) else ""


def detect_one(probe):
    """探测单项，返回结果 dict（status: ok/missing）。"""
    import glob
    paths = []
    if probe.get("vendored"):
        p = _vendored_path(probe["name"])
        paths = [p] if p else []
    else:
        for pat in probe["patterns"]:
            for root in _ROOTS:
                paths += glob.glob(os.path.join(root, pat))
        w = shutil.which(probe.get("which", "")) if probe.get("which") else ""
        if w:
            paths.append(w)
    # multi：保留全部版本；否则取最新一个
    if not paths:
        return dict(probe, path="", version="", status="missing",
                    note=DEGRADE.get(probe["name"], ""))
    if probe.get("multi"):
        seen, uniq = set(), []
        for p in paths:
            if p.lower() not in seen:
                seen.add(p.lower())
                uniq.append(p)
        paths = sorted(uniq)
    else:
        paths = [sorted(paths)[-1]]
    return dict(probe, path="; ".join(paths), version=_ver(paths[-1]),
                status="ok", note="")


def _ver(path):
    """廉价版本探测（只对极少数工具执行 --version，其余从路径/文件名猜）。"""
    if not path:
        return ""
    base = os.path.basename(os.path.dirname(os.path.dirname(path))) if "MaxQuant" in path else ""
    if "MaxQuant" in path:
        m = re.search(r"(2\.\d+\.\d+(?:\.\d+)?)", path)
        return m.group(1) if m else ""
    m = re.search(r"(\d+\.\d+(?:\.\d+)?)", path)
    return m.group(1) if m else ""


def detect_all():
    return [detect_one(p) for p in PROBES]


FIELDS = ["name", "category", "need", "status", "version", "path", "note"]


def write_report(rows, out_tsv):
    import csv
    with open(out_tsv, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def summary(rows):
    n_ok = sum(1 for r in rows if r["status"] == "ok")
    req_missing = [r["name"] for r in rows if r["status"] == "missing" and r["need"] == "required"]
    rec_missing = [r["name"] for r in rows if r["status"] == "missing" and r["need"] == "recommended"]
    return n_ok, len(rows), req_missing, rec_missing
