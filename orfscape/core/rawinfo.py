# -*- coding: utf-8 -*-
"""raw 信息提取：rawmethod（内嵌方法参数）+ rawdump（扫描表统计）+ sha256 + 汇总表。

调用约定（与 eval_ms_pipeline 一致）：
  rawmethod.exe <数据目录> <输出tsv>          —— 目录级批量，输出每 raw 一行方法摘要
  rawdump.exe  <数据目录> --all               —— 目录级批量，raw_scan_agc.tsv 写入数据目录

为不污染源数据目录，先把各分组 raw 以硬链接（同卷免管理员；失败回退拷贝）
staging 到项目临时目录再调用，产物移入项目目录后清理 staging。

子进程安全：固定参数列表 + shell=False；可执行文件为模块常量并在导入时校验存在。
"""
import csv
import hashlib
import os
import re
import shutil
import statistics
import subprocess
from pathlib import Path

from ..config import RAWDUMP_EXE, RAWMETHOD_EXE

if not (RAWMETHOD_EXE.is_file() and RAWDUMP_EXE.is_file()):
    raise FileNotFoundError("orfscape/tools 下缺少 rawmethod.exe 或 rawdump.exe，请检查安装")

# 480 模板方法段（eval_ms_pipeline 同款正则）
SEG_PAT = re.compile(
    r"CV(-?\d+):res(\d+)/(\d+):MS1_IT(\d+):MS2_IT(\d+):cyc([\d.]+)s:nce([\d.]+)")
AGC_PAT = re.compile(r"Normalised AGC Target \(%\)= (\d+)")

csv.field_size_limit(2 ** 31 - 1)


# ---------------------------------------------------------------- staging
def stage_group(raw_paths, staging_dir):
    """把一组 raw 硬链接到 staging 目录（同卷；跨卷/失败回退拷贝）。返回 staging Path。"""
    staging = Path(staging_dir)
    staging.mkdir(parents=True, exist_ok=True)
    for p in raw_paths:
        dst = staging / Path(p).name
        if dst.exists():
            continue
        try:
            os.link(p, dst)
        except OSError:
            shutil.copy2(p, dst)
    return staging


def clean_stage(staging_dir):
    shutil.rmtree(staging_dir, ignore_errors=True)


# ---------------------------------------------------------------- 工具调用
def run_rawmethod(datadir, out_tsv):
    """调 rawmethod（固定参数列表，shell=False）；返回 (ok, 输出tsv路径)。"""
    rc = subprocess.run([str(RAWMETHOD_EXE), str(datadir), str(out_tsv)],
                        shell=False, capture_output=True).returncode
    ok = rc == 0 and Path(out_tsv).exists()
    return ok, out_tsv


def run_rawdump(datadir):
    """调 rawdump --all（固定参数列表，shell=False）；返回 (ok, tsv路径)。"""
    rc = subprocess.run([str(RAWDUMP_EXE), str(datadir), "--all"],
                        shell=False, capture_output=True).returncode
    tsv = Path(datadir) / "raw_scan_agc.tsv"
    return (rc == 0 and tsv.exists()), tsv


# ---------------------------------------------------------------- 解析
def parse_method_tsv(tsv):
    """rawmethod 输出 → {raw文件名: {方法参数字段}}。

    字段：segments 列表（每段 dict: cv,res1,res2,ms1_it,ms2_it,cycle,nce）+
    汇总 resolutions/cycle_times/faims_cvs/nce/agc列表/method_len。
    """
    out = {}
    with open(tsv, newline="", encoding="utf-8", errors="replace") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            fname = (row.get("file") or "").strip()
            if not fname:
                continue
            segs = []
            for m in SEG_PAT.finditer(row.get("segments") or ""):
                segs.append({"cv": int(m.group(1)), "res1": int(m.group(2)),
                             "res2": int(m.group(3)), "ms1_it": int(m.group(4)),
                             "ms2_it": int(m.group(5)), "cycle": float(m.group(6)),
                             "nce": float(m.group(7))})
            text = (row.get("method_text") or "").replace("\\n", "\n")
            agcs = AGC_PAT.findall(text)
            out[fname] = {
                "segments": segs,
                "resolutions": row.get("resolutions", ""),
                "cycle_times": row.get("cycle_times", ""),
                "faims_cvs": row.get("faims_cvs", ""),
                "nce": row.get("nce", ""),
                "agc_pcts": agcs,
                "method_len": len(text),
            }
    return out


def parse_scan_tsv(tsv):
    """rawdump 输出 → {raw文件名: 扫描统计}。

    统计：n_scan/n_ms1/n_ms2、rt 范围、cv 集合、MS2 注入时间中位、MS1/MS2 AGC Fill 中位。
    """
    per = {}
    with open(tsv, newline="", encoding="utf-8", errors="replace") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            fname = (row.get("file") or "").strip()
            d = per.setdefault(fname, {"it_ms2": [], "fill2": [], "fill1": [],
                                       "rt": [], "cv": set(), "n_ms1": 0, "n_ms2": 0})
            try:
                order = int(row.get("msorder") or 0)
                rt = float(row.get("rt_min") or 0)
                it = float(row.get("it_ms") or 0)
                fill = float(row.get("agc_fill") or 0)
                cv = (row.get("cv") or "").strip()
            except ValueError:
                continue
            d["rt"].append(rt)
            if cv:
                d["cv"].add(cv)
            if order == 1:
                d["n_ms1"] += 1
                d["fill1"].append(fill)
            elif order == 2:
                d["n_ms2"] += 1
                d["it_ms2"].append(it)
                d["fill2"].append(fill)
    out = {}
    for fname, d in per.items():
        out[fname] = {
            "n_scan": d["n_ms1"] + d["n_ms2"],
            "n_ms1": d["n_ms1"], "n_ms2": d["n_ms2"],
            "rt_first": min(d["rt"]) if d["rt"] else "",
            "rt_last": max(d["rt"]) if d["rt"] else "",
            "cvs": ";".join(sorted(d["cv"])),
            "med_it_ms2": round(statistics.median(d["it_ms2"]), 2) if d["it_ms2"] else "",
            "med_fill_ms1": round(statistics.median(d["fill1"]), 4) if d["fill1"] else "",
            "med_fill_ms2": round(statistics.median(d["fill2"]), 4) if d["fill2"] else "",
            "_it_ms2_list": d["it_ms2"],  # 触顶率计算用（汇总后剔除）
        }
    return out


def saturation_pct(scan_stats, method_info):
    """MS2 注入时间触顶率 = count(it >= ms2_mit-0.5)/n_ms2（方法缺失返回空）。"""
    segs = (method_info or {}).get("segments") or []
    if not segs:
        return ""
    mit = max(s["ms2_it"] for s in segs) - 0.5
    its = scan_stats.get("_it_ms2_list") or []
    if not its:
        return ""
    return round(100.0 * sum(1 for x in its if x >= mit) / len(its), 2)


# ---------------------------------------------------------------- 文件信息
def sha256sum(path, chunk=8 << 20):
    """文件完整性校验和（SHA-256）。"""
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while True:
            b = fh.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def file_info(path):
    st = os.stat(path)
    import datetime
    return {"size_mb": round(st.st_size / 1048576, 1),
            "mtime": datetime.datetime.fromtimestamp(st.st_mtime).strftime("%Y-%m-%d %H:%M")}
