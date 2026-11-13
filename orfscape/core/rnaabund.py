# -*- coding: utf-8 -*-
"""转录本丰度（RNA-seq）→ 按 ORF 位点汇总。

输入是 StringTie 的**单样本** GTF（含 `cov` / `FPKM` / `TPM` 属性）；
`stringtie --merge` 的产出没有丰度，不能直接用。

口径
----
* 每个样本内：取与 ORF **同链**、且覆盖 ORF 长度 ≥ `min_frac` 的转录本，
  用其中 TPM 最大的一条代表该位点在本样本的丰度；
* 跨样本：对"有转录本覆盖"的样本取平均（下游可看 `n_sample` 判断可信度）；
* 只有 ≥1 个样本覆盖时才有值，否则记 NaN。
"""
import glob
import os
import re
from collections import defaultdict
from pathlib import Path

_ATTR = re.compile(r'(\S+)\s+"([^"]*)"')


def parse_gtf_line_attrs(field):
    return {k: v for k, v in _ATTR.findall(field)}


def read_sample_gtf(path):
    """→ [(tx, ref, contig, strand, start, end, tpm, cov)]"""
    out = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9 or c[2] != "transcript":
                continue
            at = parse_gtf_line_attrs(c[8])
            try:
                tpm = float(at.get("TPM") or "nan")
            except ValueError:
                tpm = float("nan")
            try:
                cov = float(at.get("cov") or "0")
            except ValueError:
                cov = 0.0
            if tpm != tpm:          # NaN
                continue
            tx = at.get("transcript_id") or ""
            ref = at.get("reference_id") or at.get("ref_gene_id") or tx
            try:
                s, e = int(c[3]), int(c[4])
            except ValueError:
                continue
            out.append((tx, ref, c[0], c[6], s, e, tpm, cov))
    return out


def load_samples(paths):
    """paths 可以是文件、目录或 glob。"""
    files = []
    for p in paths:
        if os.path.isdir(p):
            files.extend(sorted(glob.glob(os.path.join(p, "*.gtf"))))
        elif any(ch in p for ch in "*?["):
            files.extend(sorted(glob.glob(p)))
        else:
            files.append(p)
    files = [f for f in files if os.path.basename(f) != "merged.gtf"]
    return [(os.path.basename(f), read_sample_gtf(f)) for f in files]


def abundance_for_span(samples, contig, strand, s, e, min_frac=0.5,
                       max_len_mult=4.0, max_len=8000):
    """ORF 位点的 RNA-seq 丰度（跨样本平均）。

    为什么要长度上限：StringTie 在原核常把相邻基因连成几十 kb 的"转录本"
    （实测 chr 上 233,103–283,593 这种 50 kb 条目 TPM 上万），直接取最大 TPM
    会被这类条目主导。因此：
      * `cov`  ：按碱基加权估计 ORF 上的平均深度 Σ(cov_t × 重叠长)/ORF长（对长条目稳健）
      * `tpm`  ：仅取"基因尺度"转录本（长度 ≤ max(4×ORF, 8 kb) 且覆盖 ≥min_frac）的最大 TPM
    反义事件（同链无转录本）回退到**对侧链**——同一条 RNA 分子反向翻译即可，
    并在 `strand_used` 里标 "opposite"。

    返回 dict(cov, tpm, n_sample, tx, strand_used)
    """
    orf_len = e - s + 1
    need = max(1, int(orf_len * min_frac))
    covs, tpms, n_hit, txs = [], [], 0, set()
    used = None
    for use_strand, tag in ((strand, "same"), ("-" if strand == "+" else "+", "opposite")):
        covs, tpms, n_hit, txs = [], [], 0, set()
        for _name, recs in samples:
            depth = 0.0
            best_tpm, best_tx = 0.0, None
            for _tx, ref, ctg, st, ts, te, tpm, cov in recs:
                if ctg != contig or st != use_strand:
                    continue
                ov = min(e, te) - max(s, ts) + 1
                if ov <= 0:
                    continue
                depth += cov * ov / orf_len
                tlen = te - ts + 1
                if ov >= need and tlen <= min(max_len, max_len_mult * orf_len):
                    if tpm > best_tpm:
                        best_tpm, best_tx = tpm, ref
            if depth > 0 or best_tx:
                n_hit += 1
                covs.append(depth)
                tpms.append(best_tpm)
                if best_tx:
                    txs.add(best_tx)
        if n_hit:
            used = tag
            break
    if not n_hit:
        return dict(cov=float("nan"), tpm=float("nan"), n_sample=0, tx=[], strand_used=None)
    return dict(cov=sum(covs) / len(covs),
                tpm=(sum(tpms) / len(tpms)) if any(tpms) else float("nan"),
                n_sample=n_hit, tx=sorted(txs), strand_used=used)


# ---------------------------------------------------------------- 转录本结构 + 丰度
_ATTR2 = re.compile(r'(\S+)\s+"([^"]*)"')


def load_transcripts(merged_gtf, samples=None):
    """StringTie `--merge` 的 GTF → 转录本结构；再用单样本 GTF 的平均 TPM 标注。

    返回 [dict(id, ref, contig, strand, s, e, exons=[(s,e)...], tpm)]。
    merge 后的 GTF 没有丰度，TPM 必须从单样本 GTF 按 `reference_id` 对上。
    """
    path = os.path.realpath(merged_gtf)
    if not os.path.isfile(path):
        return []
    tx = {}
    for line in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        if not line or line.startswith("#"):
            continue
        c = line.split("	")
        if len(c) < 9 or c[2] not in ("transcript", "exon"):
            continue
        at = dict(_ATTR2.findall(c[8]))
        tid = at.get("transcript_id")
        if not tid:
            continue
        try:
            s0, e0 = int(c[3]), int(c[4])
        except ValueError:
            continue
        d = tx.get(tid)
        if d is None:
            d = tx[tid] = dict(id=tid, ref=at.get("reference_id") or at.get("ref_gene_id") or "",
                               contig=c[0], strand=c[6], s=s0, e=e0, exons=[], tpm=0.0)
        if c[2] == "exon":
            d["exons"].append((s0, e0))
        d["s"] = min(d["s"], s0)
        d["e"] = max(d["e"], e0)

    tpm_acc = defaultdict(list)
    for _name, recs in (samples or []):
        for _tx, ref, _ctg, _st, _s, _e, tpm, _cov in recs:
            if ref:
                tpm_acc[ref].append(tpm)
    for d in tx.values():
        vals = tpm_acc.get(d["id"]) or tpm_acc.get(d["ref"]) or []
        d["tpm"] = (sum(vals) / len(vals)) if vals else 0.0
        d["exons"].sort()
    return sorted(tx.values(), key=lambda d: (d["contig"], d["s"]))


# ---------------------------------------------------------------- 单样本组装模型 + 覆盖度
def assembled_models(paths, log=None, max_len=8000):
    """单样本 StringTie 组装 → 去重成转录本模型（同 contig/链/坐标合并，取最大 TPM）。

    与 `merged.gtf` 的区别：merge 会把相邻基因并成几十 kb 的条目，单样本模型才是
    "每个样本真正拼出来的转录本"。每个模型保留外显子块与平均深度 cov。
    """
    files = []
    for p in paths:
        if os.path.isdir(p):
            files.extend(sorted(glob.glob(os.path.join(p, "*.gtf"))))
        elif any(ch in p for ch in "*?["):
            files.extend(sorted(glob.glob(p)))
        else:
            files.append(p)
    files = [f for f in files if os.path.basename(f) != "merged.gtf"]

    best = {}
    dropped = 0
    for f in files:
        cur = {}
        for line in Path(f).read_text(encoding="utf-8", errors="replace").splitlines():
            if not line or line.startswith("#"):
                continue
            c = line.split("	")
            if len(c) < 9 or c[2] not in ("transcript", "exon"):
                continue
            at = dict(_ATTR2.findall(c[8]))
            tid = at.get("transcript_id")
            if not tid:
                continue
            try:
                s0, e0 = int(c[3]), int(c[4])
            except ValueError:
                continue
            d = cur.get(tid)
            if d is None:
                d = cur[tid] = dict(id=at.get("reference_id") or tid, tx=tid,
                                    contig=c[0], strand=c[6], s=s0, e=e0,
                                    exons=[], tpm=0.0, cov=0.0)
            if c[2] == "exon":
                d["exons"].append((s0, e0))
            d["s"], d["e"] = min(d["s"], s0), max(d["e"], e0)
            try:
                d["tpm"] = max(d["tpm"], float(at.get("TPM") or 0))
                d["cov"] = max(d["cov"], float(at.get("cov") or 0))
            except ValueError:
                pass
        for d in cur.values():
            if max_len and (d["e"] - d["s"] + 1) > max_len:
                # 原核基因组里 >8 kb 的"转录本"几乎都是 StringTie 把相邻基因串起来的
                dropped += 1
                continue
            d["exons"] = tuple(sorted(d["exons"])) or ((d["s"], d["e"]),)
            key = (d["contig"], d["strand"], d["exons"])
            old = best.get(key)
            if old is None or d["tpm"] > old["tpm"]:
                d["n_sample"] = 1
                best[key] = d
            else:
                old["n_sample"] += 1
    out = sorted(best.values(), key=lambda d: (d["contig"], d["s"]))
    if log:
        log("组装模型 %d 条（来自 %d 个样本 GTF；丢弃 >%d bp 的长条目 %d 条）"
            % (len(out), len(files), max_len, dropped))
    return out


def coverage_profile(models, contig, x0, x1):
    """窗口内逐碱基深度（由各模型的 cov × 覆盖区间累加，step 函数）。

    这是 StringTie 模型深度，不是 reads 的逐碱基堆叠（没有 BAM 时用它最接近）。
    """
    lo, hi = int(x0), int(x1)
    if hi <= lo:
        return [lo, hi], [0.0, 0.0]
    diff = [0.0] * (hi - lo + 2)
    for m in models:
        if m["contig"] != contig or not m.get("cov"):
            continue
        for a, b in m["exons"]:
            a2, b2 = max(a, lo), min(b, hi)
            if b2 < a2:
                continue
            diff[a2 - lo] += m["cov"]
            diff[b2 - lo + 1] -= m["cov"]
    out = []
    cur = 0.0
    for i in range(hi - lo + 1):
        cur += diff[i]
        out.append(cur)
    return list(range(lo, hi + 1)), out


def tx_support(models, contig, strand, s, e, min_frac=0.5):
    """某区间是否有组装转录本支持。

    返回 dict(n, samples, max_tpm, covered)，`covered` 表示至少一个模型覆盖
    ≥min_frac 的区间长度（这是"转录本支持"的判据，不看总量）。
    """
    need = max(1, int((e - s + 1) * min_frac))
    hit = [m for m in (models or [])
           if m["contig"] == contig and m["strand"] == strand
           and min(e, m["e"]) - max(s, m["s"]) + 1 >= need]
    if not hit:
        return dict(n=0, samples=0, max_tpm=0.0, covered=False)
    return dict(n=len(hit),
                samples=max((m.get("n_sample") or 0) for m in hit),
                max_tpm=max((m.get("tpm") or 0.0) for m in hit),
                covered=True)
