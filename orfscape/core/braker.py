# -*- coding: utf-8 -*-
"""BRAKER3 真核基因模型校正（B/C 线用）。

BRAKER3 = GeneMark-EP+（蛋白证据）+ GeneMark-ET（RNA-seq 剪接位点）+ AUGUSTUS，
输出真核基因模型（GTF + 蛋白 FASTA）。它解决的是**真核**线的问题：六框翻译
对真核基因组没有内含子概念，参考注释又可能缺基因/边界不准；BRAKER3 用 RNA-seq
与跨物种蛋白把模型校正好，再进蛋白质基因组学比对。

**它只在 Linux 跑**（GeneMark 二进制 + AUGUSTUS + Perl/BioPerl）。本模块：
  1. 生成可直接执行的运行脚本（参数全部显式，路径由调用方给）；
  2. 在能跑的环境里调用并收集日志；
  3. 解析产出 GTF → 蛋白 FASTA（供下游 L1 等价层使用）。

本机（Windows）只做前两步的准备与第 3 步的解析，运行会明确报"不可用"，
不静默跳过——避免又出现"0 结果但流程继续"。
"""
import os
import re
import subprocess

CODON = {}
_B = "TCAG"
_AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
_i = 0
for _a in _B:
    for _b in _B:
        for _c in _B:
            CODON[_a + _b + _c] = _AA[_i]
            _i += 1


def revcomp(nt):
    return nt.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def translate(nt):
    nt = nt.upper()
    return "".join(CODON.get(nt[i:i + 3], "X") for i in range(0, len(nt) - 2, 3))


def check_env(braker_cmd="braker3"):
    """BRAKER3 是否可调用。返回 (ok, 说明)。"""
    if os.name == "nt":
        return False, ("BRAKER3 需要 Linux（GeneMark + AUGUSTUS + Perl）；"
                       "本机为 Windows。请在 WSL/VM 里跑，或用 --braker-cmd 指向"
                       "可执行路径后由该环境调用。")
    from shutil import which
    p = which(braker_cmd) or (braker_cmd if os.path.isfile(braker_cmd) else None)
    if not p:
        return False, "未找到 %s（conda install -c bioconda braker3）" % braker_cmd
    return True, p


def build_run_script(path, genome, bam, protein, out_dir, threads=8,
                     species="arabidopsis", braker_cmd="braker3"):
    """写一个可直接 `bash` 的运行脚本（调用方给绝对路径）。

    BRAKER3 的输入约定：
      --genome  基因组 FASTA
      --bam     比对到该基因组的 RNA-seq BAM（可给多个，空格分隔）
      --prot_seq 跨物种蛋白 FASTA（可选，给了就走 EP+ 分支）
    """
    lines = [
        "#!/usr/bin/env bash",
        "set -euo pipefail",
        'mkdir -p "%s"' % out_dir,
        'cd "%s"' % out_dir,
        '%s \\' % braker_cmd,
        '  --genome="%s" \\' % genome,
        '  --bam="%s" \\' % bam,
    ]
    if protein:
        lines.append('  --prot_seq="%s" \\' % protein)
    lines += [
        '  --threads=%d \\' % threads,
        '  --species="%s" \\' % species,
        '  --gff3',
    ]
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(lines) + "\n")
    return path


def run(script_path, log=print, timeout=86400):
    """执行脚本。返回 (状态, 返回码)。非 Linux 环境直接返回 NOT_LINUX。"""
    if os.name == "nt":
        log("  BRAKER3 跳过：Windows 本机不可运行（请在 WSL/VM 执行生成的脚本）")
        return "NOT_LINUX", -1
    try:
        r = subprocess.run(["bash", str(script_path)], capture_output=True,
                           timeout=timeout)
    except (OSError, subprocess.SubprocessError) as exc:
        log("  BRAKER3 运行失败：%s" % exc)
        return "FAIL", -1
    return ("OK" if r.returncode == 0 else "FAIL"), r.returncode


def parse_gtf(path):
    """GTF → CDS 段列表：[(transcript_id, seqid, start1, end1, strand)]（多段=多外显子）。"""
    segs = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[2] != "CDS":
                continue
            tid = ""
            for kv in p[8].split(";"):
                kv = kv.strip()
                if kv.startswith("transcript_id") or kv.startswith("Parent"):
                    tid = kv.split()[-1].strip('"')
                    break
            segs.append((tid, p[0], int(p[3]), int(p[4]), p[6]))
    return segs


def gtf_to_protein(gtf_path, genome_path):
    """把 BRAKER GTF 的 CDS 拼接翻译成蛋白（按链方向，自动去尾 stop）。

    返回 [(protein_id, seq)]，protein_id = <tid>|<seqid>|<start>-<end>|<strand>。
    """
    genome = {}
    name = None
    buf = []
    with open(genome_path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    genome[name] = "".join(buf)
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if name:
        genome[name] = "".join(buf)

    by_tid = {}
    for tid, seqid, s, e, strand in parse_gtf(gtf_path):
        by_tid.setdefault((tid, seqid, strand), []).append((s, e))

    out = []
    for (tid, seqid, strand), segs in by_tid.items():
        nt = genome.get(seqid, "")
        if not nt:
            continue
        segs.sort()
        seq = "".join(nt[s - 1:e] for s, e in segs)
        if strand == "-":
            seq = revcomp(seq)
        aa = translate(seq)
        aa = aa.rstrip("*")
        if not aa:
            continue
        pid = "%s|%s|%d-%d|%s" % (tid, seqid, segs[0][0], segs[-1][1], strand)
        out.append((pid, aa))
    return out


def write_faa(path, recs):
    """写蛋白 FASTA（60 列折行）。path 由调用方给绝对路径。"""
    with open(path, "w", encoding="utf-8") as fh:
        for pid, seq in recs:
            fh.write(">%s\n" % pid)
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")
    return path


def filter_models(proteins, cds_ranges, min_aa=30, max_overlap=0.8):
    """去掉过短、以及与已有注释同链高度重叠的模型（与 L1 的过滤口径一致）。"""
    kept = []
    for pid, seq in proteins:
        if len(seq) < min_aa:
            continue
        m = re.search(r"\|(\d+)-(\d+)\|([+-])$", pid)
        if not m or not cds_ranges:
            kept.append((pid, seq))
            continue
        s, e, st = int(m.group(1)), int(m.group(2)), m.group(3)
        length = e - s + 1
        clash = False
        for a, b, cst in cds_ranges:
            if cst != st:
                continue
            lo, hi = max(s, a), min(e, b)
            if hi - lo + 1 > max_overlap * length:
                clash = True
                break
        if not clash:
            kept.append((pid, seq))
    return kept
