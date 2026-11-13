# -*- coding: utf-8 -*-
"""E1 已注释库搜库调度（设计书 v1.2）：pFind 常规臂 × db_annot × 每分 CV MGF。

模板与修饰库已 vendor 至 orfscape/tools/（自包含）。
cfg 生成：模板行级替换（GBK、正斜杠）；酶按组：TL=Trypsin KR _ C，
LL=LysargiNase KR _ N（pFind enzyme.ini 实名）。
断点续跑：任务目录已有 pFind-Filtered.spectra 即跳过。
子进程安全：参数列表内联 + shell=False + 在 pFind bin 目录运行（交接口径）。
"""
import os
import subprocess

_TOOLS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "tools")
_TPL_FILE = "restrict_pFind.cfg"
_MOD_FILE = "pFind_modification.ini"
ENZYME_LINE = {"TL": "enzyme=Trypsin KR _ C",
               "LL": "enzyme=LysargiNase KR _ N"}


def _local_tools_path(fname):
    p = os.path.join(_TOOLS, fname)
    if not os.path.isfile(p):
        raise FileNotFoundError(p)
    return p


def build_cfg_lines(fasta_fwd, mgf_fwd, out_dir_fwd, out_name, enzyme_code, thread=8):
    """生成 cfg 文本（字符串）。调用方负责落盘（GBK）。入参须为正斜杠绝对路径。"""
    tpl = open(_local_tools_path(_TPL_FILE), encoding="gbk").read().splitlines()
    modini = _local_tools_path(_MOD_FILE).replace("\\", "/")
    repl = {
        "thread=": "thread=" + str(int(thread)),
        "fixmod=": "fixmod=Carbamidomethyl[C];",
        "selectmod=": "selectmod=Oxidation[M];Acetyl[AnyN-term];",
        "max_clv_sites=": "max_clv_sites=2",
        "len_lower=": "len_lower=7",
        "modpath=": "modpath=" + modini,
        "fastapath=": "fastapath=" + str(fasta_fwd),
        "outputpath=": "outputpath=" + str(out_dir_fwd),
        "outputname=": "outputname=" + str(out_name),
        "msmsnum=": "msmsnum=1",
        "msmspath1=": "msmspath1=" + str(mgf_fwd),
        "enzyme=": ENZYME_LINE.get(str(enzyme_code), ENZYME_LINE["TL"]),
    }
    out_lines, seen_msms1 = [], False
    for ln in tpl:
        hit = None
        for k, v in repl.items():
            if ln.startswith(k):
                hit = (k, v)
                break
        if hit:
            out_lines.append(hit[1])
            if hit[0] == "msmspath1=":
                seen_msms1 = True
        else:
            if ln.startswith("msmspath") and seen_msms1:
                continue
            out_lines.append(ln)
    return "\n".join(out_lines) + "\n"


def run_searcher(cfg_abs, pfind_bin_dir_abs):
    """在 pFind bin 目录运行 Searcher.exe（参数列表内联 + shell=False）。"""
    searcher = os.path.join(str(pfind_bin_dir_abs), "Searcher.exe")
    if not os.path.isfile(searcher):
        raise FileNotFoundError(searcher)
    r = subprocess.run([searcher, str(cfg_abs)], shell=False, capture_output=True,
                       cwd=str(pfind_bin_dir_abs))
    return r.returncode


def read_summary(task_dir, out_name):
    """读任务目录 pFind.summary → dict（ID 率等，容错字段名）。"""
    p = os.path.join(str(task_dir), "pFind.summary")
    if not os.path.exists(p):
        alt = os.path.join(os.path.join(str(task_dir), str(out_name)), "pFind.summary")
        p = alt if os.path.exists(alt) else p
    out = {}
    if not os.path.exists(p):
        return out
    for line in open(p, encoding="gbk", errors="replace"):
        for k in ("Total MS2", "Identified MS2", "ID Rate", "Peptides", "Proteins"):
            if line.strip().startswith(k):
                out[k.replace(" ", "_")] = line.split(":", 1)[-1].strip() if ":" in line else ""
                break
    return out


def count_result_rows(spectra_file):
    """pFind-Filtered.spectra 行数（不含表头）；不存在返回 0。"""
    if not os.path.exists(str(spectra_file)):
        return 0
    with open(spectra_file, encoding="utf-8", errors="replace") as fh:
        return sum(1 for _ in fh) - 1
