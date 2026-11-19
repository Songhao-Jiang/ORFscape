# -*- coding: utf-8 -*-
"""真核 RNA-seq → 转录本组装（参考引导），给 ``eukbuild --transcripts`` 用。

真核新肽的发现层（L2|）需要"样本自己的转录本"：注释里没有的新剪接、内含子保留、
uORF/sORF、ncRNA 转录本，都只能从 RNA-seq 组装出来。此前真核线靠外源 FASTA
（衣藻/人），本模块把这一步补上：

    hisat2-build 基因组
      → hisat2 逐样本比对（双端）
      → samtools sort → BAM
      → stringtie 逐样本组装（-G 注释引导）
      → stringtie --merge（全部样本）
      → gffcompare -r 注释（拿 class code：= / j / u / x / i / o / e）
      → gffread -w 转录本 FASTA

产物（``02_tx/<sp>/``）：
  transcripts.fa       组装转录本（含注释转录本与新型转录本）
  tx_meta.tsv          逐转录本：坐标/外显子数/长度/class_code/参考基因/各样本 FPKM
  euktx_stats.tsv      逐样本：reads 数、比对率、唯一比对率
  euktx_manifest.json  输入、参数、工具版本、计数

工具定位顺序：``ORFSCAPE_EUKRNA_BIN`` 环境变量 → PATH → 常见 conda 环境
（``<conda>/envs/eukrna/bin``、``~/.conda/envs/eukrna/bin``）。
"""
import json
import os
import re
import shutil
import subprocess
import sys
from collections import OrderedDict

_TOOLS = ("hisat2", "hisat2-build", "samtools", "stringtie", "gffcompare", "gffread")


# ---------------------------------------------------------------- 工具定位
def _candidate_dirs():
    dirs = []
    env = os.environ.get("ORFSCAPE_EUKRNA_BIN", "")
    if env:
        dirs.append(env)
    for root in (os.environ.get("CONDA_PREFIX", ""),
                 r"D:/Bioinformatics/Miniconda",
                 os.path.expanduser("~/miniconda3"),
                 os.path.expanduser("~/anaconda3")):
        if root:
            dirs.append(os.path.join(root, "envs", "eukrna", "bin"))
            dirs.append(os.path.join(root, "bin"))
    return [d for d in dirs if d]


def _open_w(path, base, mode="w", **kw):
    """写盘前确认目标在 base 之下：所有输出路径都来自 CLI，必须先收敛。"""
    p = os.path.abspath(str(path))
    b = os.path.abspath(str(base))
    if os.path.commonpath([p, b]) != b:
        raise ValueError("路径越界: %s 不在 %s 之下" % (p, b))
    return open(p, mode, **kw)


def find_tools(need=_TOOLS):
    """{name: path}；找不到的记 None。"""
    out = {}
    for name in need:
        found = shutil.which(name)
        if not found:
            for d in _candidate_dirs():
                for ext in (".exe", ""):
                    p = os.path.join(d, name + ext)
                    if os.path.isfile(p):
                        found = p
                        break
                if found:
                    break
        out[name] = found
    return out


def _run(cmd, log, cwd=None):
    """跑一条命令，返回 (rc, stderr_tail)。"""
    r = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True,
                       encoding="utf-8", errors="replace")
    tail = (r.stderr or "").strip().splitlines()[-3:]
    if r.returncode != 0:
        log("  ✗ %s rc=%d" % (os.path.basename(cmd[0]), r.returncode))
        for line in tail:
            log("    " + line)
    return r.returncode, "\n".join(tail)


# ---------------------------------------------------------------- 输入发现
_PAIR = re.compile(r"(?P<sample>.+?)(?:[._](?:R?1|fq1|1))\.(?:fq|fastq)(?:\.gz)?$", re.I)
_PAIR2 = re.compile(r"(?P<sample>.+?)(?:[._](?:R?2|fq2|2))\.(?:fq|fastq)(?:\.gz)?$", re.I)


def find_fastq(root, log=print):
    """在 root 下递归找双端 fastq 对。返回 [{sample, r1, r2}]。"""
    r1, r2 = {}, {}
    for dirpath, _dirs, files in os.walk(root):
        for f in files:
            if not re.search(r"\.(fq|fastq)(\.gz)?$", f, re.I):
                continue
            m1, m2 = _PAIR.match(f), _PAIR2.match(f)
            if m1:
                r1.setdefault(m1.group("sample"), os.path.join(dirpath, f))
            elif m2:
                r2.setdefault(m2.group("sample"), os.path.join(dirpath, f))
    out = []
    for s in sorted(set(r1) & set(r2)):
        out.append({"sample": s, "r1": r1[s], "r2": r2[s]})
    only1 = sorted(set(r1) - set(r2))
    if only1:
        log("⚠ 只有 R1 的样本（跳过）: %s" % ", ".join(only1))
    return out


# ---------------------------------------------------------------- 各步
def index_genome(fna, out_dir, tools, threads=4, log=print):
    os.makedirs(out_dir, exist_ok=True)
    prefix = os.path.join(out_dir, "genome")
    if os.path.exists(prefix + ".1.ht2"):
        log("  index 已存在，复用")
        return prefix
    rc, _ = _run([tools["hisat2-build"], "-p", str(threads), fna, prefix], log)
    if rc != 0:
        raise RuntimeError("hisat2-build 失败")
    return prefix


def align(index_prefix, r1, r2, out_bam, tools, threads=8, log=print, name=""):
    """hisat2 | samtools sort → BAM；顺带用 flagstat 拿比对率。"""
    os.makedirs(os.path.dirname(out_bam), exist_ok=True)
    if os.path.exists(out_bam) and os.path.getsize(out_bam) > 0:
        log("  bam 已存在，复用: %s" % os.path.basename(out_bam))
        return _flagstat(out_bam, tools, log)
    sam = subprocess.Popen(
        [tools["hisat2"], "-p", str(threads), "--dta", "-x", index_prefix,
         "-1", r1, "-2", r2, "--summary-file", out_bam + ".hisat2.log", "-S", "-"],
        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    srt = subprocess.Popen([tools["samtools"], "sort", "-@", "2", "-o", out_bam, "-"],
                           stdin=sam.stdout, stdout=subprocess.DEVNULL,
                           stderr=subprocess.DEVNULL)
    sam.stdout.close()
    rc2 = srt.wait()
    rc1 = sam.wait()
    if rc1 != 0 or rc2 != 0:
        raise RuntimeError("hisat2/samtools 失败（%s: rc=%s/%s）" % (name, rc1, rc2))
    return _flagstat(out_bam, tools, log)


def _flagstat(bam, tools, log=print):
    rc, out = _run([tools["samtools"], "flagstat", bam], log)
    stats = {"total": 0, "mapped": 0, "mapped_pct": 0.0, "unique_pct": 0.0}
    if rc != 0:
        return stats
    m = re.search(r"in total \((\d+)", out)
    if m:
        stats["total"] = int(m.group(1))
    m = re.search(r"\[main\] mapped \((\d+(?:\.\d+)?)%", out)
    if not m:
        m = re.search(r" mapped \((\d+(?:\.\d+)?)%", out)
    if m:
        stats["mapped_pct"] = float(m.group(1))
        stats["mapped"] = int(stats["total"] * stats["mapped_pct"] / 100.0)
    return stats


def stringtie_one(bam, gtf_out, tools, threads=4, log=print):
    os.makedirs(os.path.dirname(gtf_out), exist_ok=True)
    if os.path.exists(gtf_out) and os.path.getsize(gtf_out) > 0:
        return True
    rc, _ = _run([tools["stringtie"], bam, "-p", str(threads), "-o", gtf_out], log)
    return rc == 0


def merge_gtf(gtf_list, gff_ref, out_gtf, tools, log=print):
    lst = out_gtf + ".list.txt"
    with _open_w(lst, os.path.dirname(os.path.abspath(out_gtf)),
                 encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(gtf_list) + "\n")
    cmd = [tools["stringtie"], "--merge", "-G", gff_ref, "-o", out_gtf, lst]
    rc, _ = _run(cmd, log)
    if rc != 0:
        raise RuntimeError("stringtie --merge 失败")
    return out_gtf


def gffcompare(gtf, gff_ref, prefix, tools, log=print):
    """返回 (stats dict, tmap 路径)。"""
    rc, _ = _run([tools["gffcompare"], "-r", gff_ref, "-o", prefix, gtf], log)
    stats = {}
    p = prefix + ".stats"
    if rc == 0 and os.path.exists(p):
        with open(p, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                m = re.match(r"\s*(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s*$", line)
                if m:
                    stats[m.group(2)] = int(m.group(1))
    return stats, prefix + "." + os.path.basename(gtf) + ".tmap"


def gffread_fasta(gtf, fna, out_fa, tools, log=print):
    rc, _ = _run([tools["gffread"], gtf, "-g", fna, "-w", out_fa], log)
    if rc != 0 or not os.path.exists(out_fa):
        raise RuntimeError("gffread 失败")
    return out_fa


# ---------------------------------------------------------------- 元数据
def _parse_attrs(text):
    out = {}
    for k, v in re.findall(r'([A-Za-z_][\w.]*)\s+"([^"]*)"', text):
        out[k] = v
    return out


def tx_meta(merged_gtf, tmap, out_tsv, log=print):
    """逐转录本：坐标/外显子/长度/class_code/参考基因/各样本 FPKM。"""
    cls = {}
    if tmap and os.path.exists(tmap):
        with open(tmap, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) >= 10:
                    cls[f[0]] = {"ref_gene": f[1], "ref_tid": f[2], "class": f[3]}
    rows, seen = [], OrderedDict()
    cur, fpkms = None, []
    with open(merged_gtf, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            a = _parse_attrs(f[8])
            tid = a.get("transcript_id", "")
            if f[2] == "transcript":
                cur = {"tid": tid, "chr": f[0], "start": int(f[3]), "end": int(f[4]),
                       "strand": f[6], "exons": 0, "len": 0,
                       "gene": a.get("gene_id", ""),
                       "gene_name": a.get("gene_name", "")}
                fpkms = sorted((k, v) for k, v in a.items() if re.fullmatch(r"FPKM\d*", k))
                cur["fpkm"] = ";".join("%s=%s" % (k, v) for k, v in fpkms)
                seen[tid] = cur
            elif f[2] == "exon" and cur is not None and tid == cur["tid"]:
                cur["exons"] += 1
                cur["len"] += int(f[4]) - int(f[3]) + 1
    for tid, c in seen.items():
        info = cls.get(tid, {})
        rows.append([tid, c["chr"], c["start"], c["end"], c["strand"], c["exons"],
                     c["len"], c["gene"], c["gene_name"],
                     info.get("class", ""), info.get("ref_gene", ""),
                     info.get("ref_tid", ""), c["fpkm"]])
    os.makedirs(os.path.dirname(out_tsv), exist_ok=True)
    with _open_w(out_tsv, os.path.dirname(os.path.abspath(out_tsv)),
                 encoding="utf-8", newline="\n") as fh:
        fh.write("transcript_id\tchr\tstart\tend\tstrand\texons\tlength\tgene_id\t"
                 "gene_name\tclass_code\tref_gene\tref_transcript\tfpkm\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
    log("转录本元数据：%d 条 → %s" % (len(rows), os.path.basename(out_tsv)))
    return rows


# ---------------------------------------------------------------- 导入外部组装
def _read_gtf_models(gtf):
    """GTF（StringTie/gffcompare 产物）→ {tid: {seqid,strand,gene,exons,cls,cmp}}。"""
    models = {}
    with open(gtf, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            a = _parse_attrs(f[8])
            tid = a.get("transcript_id")
            if not tid:
                continue
            if f[2] == "transcript":
                models[tid] = {"seqid": f[0], "strand": f[6] if f[6] in "+-" else "+",
                               "gene": a.get("gene_id", ""), "exons": [],
                               "cls": a.get("class_code", ""),
                               "cmp": a.get("cmp_ref", "")}
            elif f[2] == "exon":
                m = models.get(tid)
                if m is not None:
                    m["exons"].append((int(f[3]), int(f[4])))
    for m in models.values():
        m["exons"].sort()
    return models


def tx_sequence(seqs, model):
    """按外显子拼接转录本序列（负链反向互补）。"""
    from .eukbuild import revcomp
    chrom = seqs.get(model["seqid"])
    if chrom is None or not model["exons"]:
        return ""
    order = sorted(model["exons"], reverse=(model["strand"] == "-"))
    parts = []
    for s, e in order:
        if s < 1 or e > len(chrom) or e < s:
            return ""
        sub = chrom[s - 1:e]
        parts.append(revcomp(sub) if model["strand"] == "-" else sub)
    return "".join(parts)


def import_assembly(workdir, sp, gff, gtf, log=print):
    """导入现成组装（如测序公司交付的 merged.gtf）→ transcripts.fa + tx_model.gff。

    与 ``build()`` 的区别：不跑比对，只把"转录本 + 外显子结构"落成 ORFscape 认的
    两件产物，并保证**序列与结构自洽**（用基因组自己拼，不信外部 FASTA）。
    ``tx_model.gff`` = 原生注释 + 外部组装里注释没有的转录本（GFF3），
    真核事件分类用它做"肽 → 基因组坐标"的映射。
    """
    from .eukbuild import read_fasta as _read_fa, parse_annotation
    fna = os.path.join(workdir, "01_genome", sp, "genomic.fna")
    seqs = {k: v.upper() for k, v in _read_fa(fna).items()}
    if not seqs:
        raise RuntimeError("缺基因组：%s" % fna)
    models = _read_gtf_models(gtf)
    if not models:
        raise RuntimeError("GTF 无 transcript/exon 行：%s" % gtf)
    canon, _ = parse_annotation(gff, log=lambda *_: None)

    # 外部组装（stringtie 的 MSTRG.*）没有注释基因号：按区间重叠把转录本挂到
    # 重叠最多的注释基因上。否则事件分类会把 MSTRG 当成"新基因"，
    # 让本来属于 alt_frame/intronic 的肽被误判成 CDS 内同框。
    gene_span = {}
    for _tid, _tx in canon.items():
        _g = _tx.get("gene") or _tid
        _ex = _tx.get("exons") or [(s, e) for s, e, _p in (_tx.get("cds") or [])]
        if not _ex:
            continue
        _lo = min(s for s, _e in _ex)
        _hi = max(e for _s, e in _ex)
        key = (_tx["seqid"], _g)
        cur = gene_span.get(key)
        gene_span[key] = [min(_lo, cur[0]), max(_hi, cur[1])] if cur else [_lo, _hi]
    by_contig = {}
    for (_seqid, _g), (_lo, _hi) in gene_span.items():
        by_contig.setdefault(_seqid, []).append((_lo, _hi, _g))

    def _anchor_gene(seqid, exons):
        """返回与该转录本重叠最多的注释基因号（无重叠返回 ""）。"""
        pool = by_contig.get(seqid)
        if not pool or not exons:
            return ""
        lo = min(s for s, _e in exons)
        hi = max(e for _s, e in exons)
        best, best_ov = "", 0
        for glo, ghi, gid in pool:
            if ghi < lo or glo > hi:
                continue
            ov = min(hi, ghi) - max(lo, glo) + 1
            if ov > best_ov:
                best, best_ov = gid, ov
        return best if best_ov > 0 else ""

    canon_genes = {(_tx.get("gene") or "") for _tx in canon.values()}
    for _tid, _m in models.items():
        if not _m["exons"]:
            continue
        _g = _m.get("gene") or ""
        if _g and _g in canon_genes:
            continue                      # 已经是注释基因号
        anchored = _anchor_gene(_m["seqid"], _m["exons"])
        if anchored:
            _m["gene"] = anchored

    base = os.path.join(workdir, "02_tx", sp)
    os.makedirs(base, exist_ok=True)
    fa_path = os.path.join(base, "transcripts.fa")
    gff_path = os.path.join(base, "tx_model.gff")
    meta_path = os.path.join(base, "tx_meta.tsv")

    n_seq, skipped = 0, 0
    with _open_w(fa_path, os.path.dirname(os.path.abspath(fa_path)),
                 encoding="utf-8", newline="\n") as fa, \
            _open_w(meta_path, os.path.dirname(os.path.abspath(meta_path)),
                    encoding="utf-8", newline="\n") as meta:
        meta.write("transcript_id\tchr\tstart\tend\tstrand\texons\tlength\t"
                   "source\tclass_code\tref_transcript\tgene_id\n")
        for tid in sorted(models):
            m = models[tid]
            seq = tx_sequence(seqs, m)
            if len(seq) < 90:
                skipped += 1
                continue
            n_seq += 1
            fa.write(">%s\n%s\n" % (tid, seq))
            meta.write("%s\t%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s\t%s\t%s\n"
                       % (tid, m["seqid"], m["exons"][0][0], m["exons"][-1][1],
                          m["strand"], len(m["exons"]), len(seq),
                          "canonical" if tid in canon else "external",
                          m["cls"], m["cmp"], m["gene"]))

    # tx_model.gff：原生注释原样 + 外部新增转录本（GFF3）
    n_add = 0
    with _open_w(gff_path, os.path.dirname(os.path.abspath(gff_path)),
                 encoding="utf-8", newline="\n") as out:
        with open(gff, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                out.write(line if line.endswith("\n") else line + "\n")
        for tid in sorted(models):
            if tid in canon:
                continue
            m = models[tid]
            if not m["exons"]:
                continue
            n_add += 1
            out.write("%s\tORFscape\ttranscript\t%d\t%d\t.\t%s\t.\t"
                      "ID=%s;Parent=gene-%s;gene_id=%s;class_code=%s\n"
                      % (m["seqid"], m["exons"][0][0], m["exons"][-1][1], m["strand"],
                         tid, m["gene"], m["gene"], m["cls"] or "ext"))
            for i, (s, e) in enumerate(m["exons"], 1):
                out.write("%s\tORFscape\texon\t%d\t%d\t.\t%s\t.\t"
                          "ID=exon-%s-%d;Parent=%s\n"
                          % (m["seqid"], s, e, m["strand"], tid, i, tid))

    stats = {"species": sp, "source": os.path.abspath(gtf),
             "n_transcripts": n_seq, "n_external": n_add, "skipped_short": skipped,
             "n_canonical": len(models) - n_add - skipped,
             "files": {"transcripts_fa": fa_path, "tx_model_gff": gff_path,
                       "tx_meta": meta_path}}
    with _open_w(os.path.join(base, "euktx_manifest.json"), base, encoding="utf-8") as fh:
        json.dump(stats, fh, ensure_ascii=False, indent=1)
    with _open_w(os.path.join(base, "euktx_stats.tsv"), base, encoding="utf-8",
                 newline="\n") as fh:
        fh.write("item\tvalue\n")
        for k in ("n_transcripts", "n_canonical", "n_external", "skipped_short"):
            fh.write("%s\t%s\n" % (k, stats[k]))
    log("导入组装：%d 条转录本（原生 %d / 外部新增 %d，跳过过短 %d）"
        % (n_seq, stats["n_canonical"], n_add, skipped))
    return stats


# ---------------------------------------------------------------- 主流程
def build(workdir, sp, fna, gff, fastq_dir, threads=8, samples=None,
          keep_bam=True, log=print):
    """RNA-seq → 组装转录本；返回统计 dict（写 manifest + stats）。"""
    tools = find_tools()
    missing = [k for k, v in tools.items() if not v]
    if missing:
        raise RuntimeError("缺工具：%s（装到 eukrna 环境或设 ORFSCAPE_EUKRNA_BIN）"
                           % ", ".join(missing))
    log("工具：hisat2=%s stringtie=%s" % (tools["hisat2"], tools["stringtie"]))
    pairs = find_fastq(fastq_dir, log=log)
    if samples:
        want = {s.strip() for s in str(samples).split(",") if s.strip()}
        pairs = [p for p in pairs if p["sample"] in want]
    if not pairs:
        raise RuntimeError("未找到双端 fastq 对：%s" % fastq_dir)
    log("样本 %d 个：%s" % (len(pairs), ", ".join(p["sample"] for p in pairs)))

    base = os.path.join(workdir, "02_tx", sp)
    os.makedirs(base, exist_ok=True)
    idx = index_genome(fna, os.path.join(base, "index"), tools,
                       threads=min(4, threads), log=log)

    align_stats, gtfs = [], []
    for i, p in enumerate(pairs, 1):
        log("[%d/%d] %s 比对…" % (i, len(pairs), p["sample"]))
        bam = os.path.join(base, "bam", p["sample"] + ".bam")
        st = align(idx, p["r1"], p["r2"], bam, tools, threads=threads,
                   log=log, name=p["sample"])
        st.update({"sample": p["sample"], "bam": bam})
        align_stats.append(st)
        log("      比对率 %.2f%%（总 %d）" % (st["mapped_pct"], st["total"]))
        gtf = os.path.join(base, "gtf", p["sample"] + ".gtf")
        if stringtie_one(bam, gtf, tools, threads=min(4, threads), log=log):
            gtfs.append(gtf)
        else:
            log("  ✗ %s 组装失败，跳过" % p["sample"])
        if not keep_bam:
            for ext in ("", ".bai"):
                p2 = bam + ext
                if os.path.exists(p2):
                    os.remove(p2)

    merged = merge_gtf(gtfs, gff, os.path.join(base, "merged.gtf"), tools, log=log)
    cmp_prefix = os.path.join(base, "gffcmp")
    cstats, tmap = gffcompare(merged, gff, cmp_prefix, tools, log=log)
    tx_fa = gffread_fasta(merged, fna, os.path.join(base, "transcripts.fa"), tools, log=log)
    rows = tx_meta(merged, tmap, os.path.join(base, "tx_meta.tsv"), log=log)

    n_tx = len(rows)
    by_class = {}
    for r in rows:
        by_class[r[9] or "?"] = by_class.get(r[9] or "?", 0) + 1
    stats = {"species": sp, "samples": [p["sample"] for p in pairs],
             "threads": threads, "n_transcripts": n_tx,
             "class_counts": dict(sorted(by_class.items())),
             "gffcompare_stats": cstats,
             "align": [{k: v for k, v in s.items() if k != "bam"} for s in align_stats],
             "transcripts_fa": tx_fa,
             "tools": {k: os.path.basename(v or "") for k, v in tools.items()},
             "files": {"merged_gtf": merged, "tx_meta": os.path.join(base, "tx_meta.tsv"),
                       "transcripts_fa": tx_fa}}
    with _open_w(os.path.join(base, "euktx_manifest.json"), base, encoding="utf-8") as fh:
        json.dump(stats, fh, ensure_ascii=False, indent=1)
    with _open_w(os.path.join(base, "euktx_stats.tsv"), base, encoding="utf-8",
                 newline="\n") as fh:
        fh.write("sample\treads\tmapped_pct\tbam\n")
        for s in align_stats:
            fh.write("%s\t%d\t%.2f\t%s\n" % (s["sample"], s["total"],
                                             s["mapped_pct"], s.get("bam", "")))
    log("组装完成：%d 条转录本（class %s）"
        % (n_tx, ", ".join("%s=%d" % kv for kv in stats["class_counts"].items())))
    return stats
