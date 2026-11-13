# -*- coding: utf-8 -*-
"""D0 电压拆分（设计书 v1.2）：多 CV raw → 分 CV MGF（pFind dta 风格标题）+ 对账。

方案（2026-10-31 split_test 实战定稿）：
  msconvert <raw> --mgf --filter "peakPicking true" --filter "scanEvent <ev_ms2>" \
    --filter "titleMaker <stem>_<cv>.<ScanNumber>.<ScanNumber>.<ChargeState>"
事件映射（480 双 CV 模板方法）：scanEvent 2=MS2@CV-45，4=MS2@CV-65——不写死：
以 rawdump 扫描表 cv 列的 MS2 计数做对账，失配即标记 FAILED。
断点续跑：目标 MGF 存在且谱数对账通过则跳过。
子进程安全：固定参数列表 + shell=False + 可执行文件存在性校验（同 rawinfo.py 风格）。
"""
import csv
import os
import re
import subprocess
from collections import Counter


def _count_ions(mgf):
    n = 0
    with open(mgf, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("BEGIN IONS"):
                n += 1
    return n


def scan_ms2_by_cv(scan_tsv):
    """rawdump 扫描表 → {(文件名, cv字符串): MS2 计数}。cv 列如 -45.000。"""
    per = Counter()
    with open(scan_tsv, newline="", encoding="utf-8", errors="replace") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if (row.get("msorder") or "") == "2":
                cv = (row.get("cv") or "").strip()
                fname = (row.get("file") or "").strip()
                if cv and fname:
                    per[(fname, cv)] += 1
    return dict(per)


def build_rows(raw_info_tsv):
    """raw_info.tsv → 待拆分行列表（needs_cv_split=1 且非 superseded）。"""
    rows = []
    with open(raw_info_tsv, newline="", encoding="utf-8-sig") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r.get("needs_cv_split") == "1" and r.get("superseded", "") != "yes":
                rows.append(r)
    return rows


def split_one(raw_path, out_mgf, cv_tag, scan_event, msconvert_exe, stem):
    """单 raw 单 CV 拆分。返回 (状态, 谱数)。

    参数列表内联 + shell=False；title 过滤器参数为受控构造（cv_tag/scan_event
    仅接受调用方常量，stem 仅取文件名安全字符）。
    """
    if not os.path.isfile(str(msconvert_exe)):
        return "EXE_MISSING", 0
    title = f"{stem}_{cv_tag}.<ScanNumber>.<ScanNumber>.<ChargeState>"
    r = subprocess.run(
        [str(msconvert_exe), str(raw_path), "--mgf",
         "--outdir", str(os.path.dirname(str(out_mgf))),
         "--outfile", str(os.path.basename(str(out_mgf))),
         "--filter", "peakPicking true",
         "--filter", f"scanEvent {int(scan_event)}",
         "--filter", f"titleMaker {title}"],
        shell=False, capture_output=True)
    if r.returncode != 0 or not os.path.exists(out_mgf):
        return "MSCONVERT_FAIL", 0
    n = _count_ions(out_mgf)
    return ("OK", n) if n > 0 else ("EMPTY", n)


def convert_pre_split(raw_path, out_mgf, cv_tag, msconvert_exe, stem):
    """预拆分 raw（一个 raw 一个 CV）→ MGF：同 split_one，但不加 scanEvent 过滤。

    返回 (状态, 谱数)。参数列表内联 + shell=False；stem/cv_tag 仅取安全字符。
    """
    if not os.path.isfile(str(msconvert_exe)):
        return "EXE_MISSING", 0
    title = f"{stem}_{cv_tag}.<ScanNumber>.<ScanNumber>.<ChargeState>"
    r = subprocess.run(
        [str(msconvert_exe), str(raw_path), "--mgf",
         "--outdir", str(os.path.dirname(str(out_mgf))),
         "--outfile", str(os.path.basename(str(out_mgf))),
         "--filter", "peakPicking true",
         "--filter", f"titleMaker {title}"],
        shell=False, capture_output=True)
    if r.returncode != 0 or not os.path.exists(out_mgf):
        return "MSCONVERT_FAIL", 0
    n = _count_ions(out_mgf)
    return ("OK", n) if n > 0 else ("EMPTY", n)


_PPARSE_PARA = """[Basic Options]
datapath = {datapath}

[Advanced Options]
co-elute = 0
input_format = raw
isolation_width = 2
mars_threshold = -0.34
ipv_file = .\\IPV.txt
trainingset = EmptyPath

[Internal Switches]
output_mars_y = 0
delete_msn = 0
output_mgf = 1
output_pf = 1
debug_mode = 0
check_activationcenter = 1
output_all_mars_y = 0
rewrite_files = 0
export_unchecked_mono = 0
cut_similiar_mono = 1
mars_model = 4
output_trainingdata = 0
"""


def pparse_para_text(datapath):
    """pParse 参数文件正文（调用方自己写盘，datapath 传 raw 所在目录）。"""
    return _PPARSE_PARA.format(datapath=datapath)


def run_pparse(pfind_bin, para_path, mgf_expect):
    """用外部写好的 para 调 pParse，并核对产出的 MGF。

    pParse 必须从 pFind 的 bin 目录运行（同目录要有 `xtract_raw.exe` 与
    `Mono.Options.dll`），因此 `cwd=pfind_bin`。MGF 由 pParse 写在 raw 同目录，
    名字为 `<raw stem>_HCDFT.mgf`——`mgf_expect` 传这个预期路径，函数只做核对。

    与 msconvert 的差别：TITLE 末尾带 `.0.dta`（pFind 原生命名）。**同一项目内
    不要混用两种转换器**，否则 File_Name 对不上；下游 `split_title` 会剥 `.dta`。
    返回 (状态, 谱数)。
    """
    exe = os.path.join(str(pfind_bin), "pParse.exe")
    if not os.path.isfile(exe):
        return "PPARSE_MISSING", 0
    try:
        r = subprocess.run([exe, str(para_path)], cwd=str(pfind_bin),
                           shell=False, capture_output=True, timeout=3600)
    except (OSError, subprocess.SubprocessError):
        return "PPARSE_FAIL", 0
    if r.returncode != 0:
        return "PPARSE_FAIL", 0
    if not os.path.exists(str(mgf_expect)):
        return "PPARSE_FAIL", 0
    n = _count_ions(str(mgf_expect))
    return ("OK", n) if n > 0 else ("EMPTY", n)


def cv_from_tag(tag):
    """'_45' → '-45.000'（rawdump cv 键归一）。"""
    v = tag.lstrip("_")
    return f"-{v}.000" if not v.startswith("-") else f"{v}.000"


def reconcile(mgf_count, expect_map, raw_name, cv_key):
    """对账：expect_map 以 (raw文件名, cv) 为键。"""
    return expect_map.get((raw_name, cv_key)) == mgf_count


TITLE_RE = re.compile(r"^TITLE=(\S+)\.(\d+)\.\d+\.(\d+)\s*$")


def mgf_cv_stem(mgf):
    """校验 MGF 内标题 stem 是否统一（防 titleMaker 未生效）。"""
    stems = set()
    with open(mgf, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("TITLE="):
                m = TITLE_RE.match(line.strip())
                if m:
                    stems.add(m.group(1))
                if len(stems) > 1:
                    return False
    return len(stems) <= 1
