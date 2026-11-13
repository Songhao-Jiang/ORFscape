# -*- coding: utf-8 -*-
"""L6 转录本层管线（Windows 原生，零 Linux 工具依赖）。

路线：ENA 直下 FASTQ → 采样 BLASTN 比对基因组 → 覆盖分析找转录区域
     → 过滤已注释基因 → orfipy 六框 ORF → 新肽过滤 → L6_tx_orfs.faa

替代 VM 方案（hisat2/samtools/StringTie/rnaSPAdes 均不装）。
适用于原核（无剪接：转录区域=连续基因组覆盖区段）。
"""
import csv
import gzip
import io
import os
import subprocess
import urllib.request

from .. import localpaths

# ENA 域名白名单（Mimosa 约束）
_ENA = "https://ftp.sra.ebi.ac.uk"

def ena_fastq_urls(run_acc):
    """查询 ENA → [(url_R1, url_R2)]。"""
    import json
    url = (f"https://www.ebi.ac.uk/ena/portal/api/filereport?"
           f"accession={run_acc}&result=read_run&fields=fastq_ftp&format=json")
    try:
        resp = urllib.request.urlopen(url, timeout=15)
        data = json.loads(resp.read())
        ftp = data[0].get("fastq_ftp", "") if data else ""
        if ftp:
            f1, f2 = ftp.split(";")
            return ("https://" + f1, "https://" + f2)
    except Exception:
        pass
    return None


def download_fastq(run_acc, out_dir, max_mb=500):
    """ENA 直下 FASTQ（.fastq.gz），返回 [R1 路径, R2 路径]。"""
    os.makedirs(out_dir, exist_ok=True)
    urls = ena_fastq_urls(run_acc)
    if not urls:
        return None
    paths = []
    for i, url in enumerate(urls, 1):
        fname = f"{run_acc}_{i}.fastq.gz"
        fpath = os.path.join(out_dir, fname)
        if os.path.exists(fpath) and os.path.getsize(fpath) > 1000:
            paths.append(fpath)
            continue
        try:
            urllib.request.urlretrieve(url, fpath)
            paths.append(fpath)
        except Exception:
            pass
    return paths if len(paths) == 2 else paths


def sample_reads_to_fasta(fastq_gz, out_fasta, n_reads=100000):
    """从 FASTQ.gz 采样前 n 条读段写 FASTA（单端）。"""
    count = 0
    with gzip.open(fastq_gz, "rt") as fi, open(out_fasta, "w") as fo:
        for i, line in enumerate(fi):
            if i % 4 == 0 and line.startswith("@"):
                if count >= n_reads:
                    break
                fo.write(f">r{count}\n")
                count += 1
            elif i % 4 == 1:
                fo.write(line)
    return count


def blastn_reads(reads_fasta, genome_fna, out_tsv, blastn_exe, max_hsps=1):
    """BLASTN megablast 读段 vs 基因组。"""
    cmd = [str(blastn_exe), "-query", str(reads_fasta), "-subject", str(genome_fna),
           "-task", "megablast", "-outfmt", "6", "-evalue", "1e-10",
           "-max_hsps", str(max_hsps), "-num_threads", "4",
           "-out", str(out_tsv)]
    r = subprocess.run(cmd, shell=False, capture_output=True)
    return r.returncode == 0


def coverage_regions(blast_tsv, min_depth=5, min_len=60, merge_gap=100):
    """BLASTN 输出 → 覆盖 ≥min_depth 的连续区域列表 [(start, end, depth)]。"""
    from collections import defaultdict
    depth = defaultdict(int)
    with open(blast_tsv) as fh:
        for line in fh:
            p = line.split("\t")
            if len(p) >= 9:
                s, e = int(p[8]), int(p[9])
                lo, hi = min(s, e), max(s, e)
                for pos in range(lo, hi + 1):
                    depth[pos] += 1
    # 找 ≥min_depth 的连续区域
    regions = []
    in_region = False
    start = 0
    for pos in sorted(depth):
        d = depth[pos]
        if d >= min_depth and not in_region:
            start = pos
            in_region = True
        elif d < min_depth and in_region:
            if pos - start >= min_len:
                regions.append((start, pos, sum(depth[p] for p in range(start, pos))))
            in_region = False
    if in_region and max(depth) - start >= min_len:
        regions.append((start, max(depth), 0))
    # 合并近邻区域
    merged = []
    for r in regions:
        if merged and r[0] - merged[-1][1] <= merge_gap:
            merged[-1] = (merged[-1][0], r[1], merged[-1][2] + r[2])
        else:
            merged.append(r)
    return merged


def extract_genomic_regions(genome_fna, regions, out_fasta):
    """从基因组提取指定区域序列。"""
    # 读取基因组
    seqs = {}
    name, buf = "", []
    for line in open(genome_fna):
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(buf)
            name = line[1:].strip().split()[0]
            buf = []
        else:
            buf.append(line.strip().upper())
    if name:
        seqs[name] = "".join(buf)
    # 只取第一个 contig（细菌通常单染色体）
    chr_name = list(seqs.keys())[0]
    chr_seq = seqs[chr_name]
    with open(out_fasta, "w") as fo:
        for i, (s, e, d) in enumerate(regions):
            fo.write(f">txregion_{i}_{chr_name}_{s}_{e}\n")
            fo.write(chr_seq[s:e] + "\n")
    return len(regions)


def run_orfipy(input_fasta, out_pep, orfipy_exe):
    """orfipy 六框 ORF 调用。"""
    cmd = [str(orfipy_exe), str(input_fasta), "--pep", str(out_pep),
           "--min", "30", "--start-codons", "ATG,GTG,TTG,CTG"]
    r = subprocess.run(cmd, shell=False, capture_output=True)
    return r.returncode == 0


def _blastn_exe(tools_dir=""):
    """BLASTN 可执行文件：ORFSCAPE_BLASTN / paths.local.json 优先，其次 PATH。"""
    hit = localpaths.exe("blastn", r"D:\Bioinformatics\BLAST\bin", name="blastn")
    if hit:
        return hit
    cand = os.path.join(str(tools_dir or ""), "blastn.exe")
    return cand if os.path.isfile(cand) else ""


def build_l6(species_code, genome_fna, sra_tsv, work_dir, tools_dir, log=print):
    """L6 层全流程。"""
    import tempfile
    from .dbbuild import read_faa, il_norm, digest
    tx_dir = os.path.join(work_dir, "02_tx", species_code)
    os.makedirs(tx_dir, exist_ok=True)
    fastq_dir = os.path.join(tx_dir, "fastq")
    blastn = _blastn_exe(tools_dir)
    if not blastn:
        log(f"  [{species_code}] 缺 BLASTN：装 BLAST+ 或设 ORFSCAPE_BLASTN，L6 跳过")
        return None
    orfipy = os.path.join(os.path.expanduser("~"), "AppData", "Roaming",
                          "Python", "Python310", "Scripts", "orfipy.exe")
    # 1. 下载 + 采样 + 比对
    all_blast = []
    runs = []
    with open(sra_tsv, newline="", encoding="utf-8-sig") as fh:
        runs = [r["run"] for r in csv.DictReader(fh, delimiter="\t")]
    for i, run in enumerate(runs[:6], 1):  # 取前 6 个 run（控制时间）
        log(f"  [{species_code}] run {i}/{min(len(runs),6)}: {run}")
        fq = download_fastq(run, fastq_dir)
        if not fq:
            log(f"    下载失败，跳过")
            continue
        fa = os.path.join(tx_dir, f"{run}_sample.fa")
        n = sample_reads_to_fasta(fq[0], fa, n_reads=50000)
        log(f"    采样 {n} 读段")
        bt = os.path.join(tx_dir, f"{run}_blast.tsv")
        if blastn_reads(fa, genome_fna, bt, blastn):
            all_blast.append(bt)
            log(f"    BLASTN ✓")
        else:
            log(f"    BLASTN ✗")
    if not all_blast:
        log(f"  [{species_code}] 无 BLASTN 结果，L6 跳过")
        return None
    # 2. 合并覆盖区域
    all_regions = []
    for bt in all_blast:
        all_regions.extend(coverage_regions(bt, min_depth=3, min_len=90))
    # 去重合并
    all_regions.sort()
    merged = []
    for r in all_regions:
        if merged and r[0] <= merged[-1][1] + 50:
            merged[-1] = (merged[-1][0], max(merged[-1][1], r[1]), 0)
        else:
            merged.append(list(r))
    log(f"  [{species_code}] 转录覆盖区域 {len(merged)} 个")
    # 3. 提取 + ORF
    regions_fa = os.path.join(tx_dir, "transcribed_regions.fa")
    n_reg = extract_genomic_regions(genome_fna, merged, regions_fa)
    orfs_pep = os.path.join(tx_dir, "L6_tx_orfs.faa")
    if orfipy and os.path.isfile(orfipy):
        ok = run_orfipy(regions_fa, orfs_pep, orfipy)
        if ok:
            n_orf = sum(1 for l in open(orfs_pep) if l.startswith(">"))
            log(f"  [{species_code}] L6 ORF {n_orf} 条 → {orfs_pep}")
            return orfs_pep
    log(f"  [{species_code}] orfipy 不可用或失败")
    return None
