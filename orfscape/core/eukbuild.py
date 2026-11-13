# -*- coding: utf-8 -*-
"""真核建库：GFF3/GTF + 基因组 → 剪接后的蛋白层 + 剪接连接肽层。

原核那一套（六框 / 注释 CDS 直取）在真核会错，因为基因是**断裂**的：

- 注释蛋白必须把同一转录本的 CDS 按坐标顺序拼起来、负链先反向互补，再翻译；
- 新肽证据除"新基因/新读框"外，还多一类**剪接连接肽**（跨外显子边界的肽），
  必须显式造进库里，pFind 才可能命中；
- 转录本层（RNA-seq 组装）按 3 框翻译，替代原核的六框层。

本模块只做"建库"这一段，产出的 FASTA 直接可以进 E1/E2 搜索。

    orfscape eukbuild --workdir <proj> --sp SCER \
        --fna genome.fna --gff genomic.gff [--transcripts tx.fa]

产物（``02_db/<sp>/``）：
  L0_annot.faa      剪接后的注释蛋白
  L0_junction.faa   剪接连接肽（跨外显子，两侧各 N 个 aa）
  L2_txorf.faa      转录本 3 框 ORF（给了 --transcripts 才有）
  euk_stats.tsv     逐转录本统计（外显子数/长度/连接点数/丢弃原因）
  euk_manifest.json 输入与参数、计数
"""
import gzip
import json
import os
import re
from collections import defaultdict

# ---------------------------------------------------------------- 序列工具
_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq):
    return seq.translate(_COMP)[::-1]


_BASES = "TCAG"
# NCBI 密码子表的碱基顺序是 T,C,A,G；三张用得到的表：
#   1 标准；2 脊椎动物线粒体（TGA=W, ATA=M, AGA/AGG=*）；
#   3 酵母线粒体（TGA=W, ATA=M, CTN=T）
_AMINO = {
    1: "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    2: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",
    3: "FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
}


def _build_table(spec):
    return {a + b + c: spec[i * 16 + j * 4 + k]
            for i, a in enumerate(_BASES)
            for j, b in enumerate(_BASES)
            for k, c in enumerate(_BASES)}


CODONS = {tid: _build_table(spec) for tid, spec in _AMINO.items()}
CODON = CODONS[1]


def translate(nt, table=1, stop_at_stop=False):
    """按指定密码子表翻译（1 标准 / 2 脊椎线粒体 / 3 酵母线粒体）；非 ACGT 记 X。"""
    codons = CODONS.get(table) or CODONS[1]
    out = []
    n = len(nt) - len(nt) % 3
    for i in range(0, n, 3):
        aa = codons.get(nt[i:i + 3].upper(), "X")
        if aa == "*" and stop_at_stop:
            break
        out.append(aa)
    return "".join(out)


def read_fasta(path):
    """{id: sequence}；id 取第一个空白前的部分。"""
    out = {}
    name, buf = None, []
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name:
                    out[name] = "".join(buf)
                name = line[1:].split()[0]
                buf = []
            elif name:
                buf.append(line.strip())
    if name:
        out[name] = "".join(buf)
    return out


# ---------------------------------------------------------------- 注释解析
_RE_GFF = re.compile(r"([A-Za-z_][\w.]*)=([^;]*)")
_RE_GTF = re.compile(r'([A-Za-z_][\w.]*)\s+"([^"]*)"')


def _attrs(text, fmt):
    out = {}
    if fmt == "gtf":
        for k, v in _RE_GTF.findall(text):
            out[k] = v
    else:
        for k, v in _RE_GFF.findall(text):
            out[k] = v.replace("%3B", ";").replace("%2C", ",").replace("%3D", "=")
    return out


def parse_annotation(path, log=print):
    """解析 GFF3/GTF。

    返回 ``(txs, stats)``；``txs`` 形如::

        {"XM_1": {"seqid": "NC_1", "strand": "+", "gene": "G1",
                  "biotype": "mRNA", "cds": [(s, e, phase)], "exons": [(s, e)]}}

    坐标一律 1-based 闭区间（GFF 约定）。
    """
    fmt = "gtf" if str(path).lower().endswith((".gtf", ".gtf.gz")) else "gff3"
    opener = gzip.open if str(path).endswith(".gz") else open
    txs = {}
    cds = defaultdict(list)
    exons = defaultdict(list)
    tables = {}                      # 转录本 → 密码子表（线粒体用 2/3）
    stats = {"lines": 0, "cds": 0, "exon": 0, "tx": 0, "no_parent": 0}
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            stats["lines"] += 1
            seqid, _src, feat, start, end, _sc, strand, phase, attr = f[:9]
            a = _attrs(attr, fmt)
            if feat in ("mRNA", "transcript"):
                tid = a.get("ID") or a.get("transcript_id")
                if not tid:
                    continue
                parent = a.get("Parent") or a.get("gene_id") or ""
                txs[tid] = {"seqid": seqid,
                            "strand": strand if strand in "+-" else "+",
                            "gene": parent.split(",")[0],
                            "biotype": a.get("biotype") or a.get("transcript_biotype") or "",
                            "cds": [], "exons": [], "table": 1}
                stats["tx"] += 1
            elif feat in ("CDS", "exon"):
                tid = a.get("Parent") or a.get("transcript_id")
                if not tid:
                    stats["no_parent"] += 1
                    continue
                ph = 0
                try:
                    ph = int(phase)
                except (TypeError, ValueError):
                    ph = 0
                try:
                    tbl = int(a.get("transl_table") or 1)
                except (TypeError, ValueError):
                    tbl = 1
                for one in tid.split(","):
                    if feat == "CDS":
                        cds[one].append((int(start), int(end), ph))
                        if tbl != 1:
                            tables[one] = tbl
                        stats["cds"] += 1
                    else:
                        exons[one].append((int(start), int(end)))
                        stats["exon"] += 1
    for tid, tbl in tables.items():
        if tid in txs:
            txs[tid]["table"] = tbl
    for tid, segs in cds.items():
        if tid not in txs:                       # CDS 有、mRNA 行缺失时兜底
            txs[tid] = {"seqid": "", "strand": "+", "gene": tid, "biotype": "",
                        "cds": [], "exons": [], "table": tables.get(tid, 1)}
        txs[tid]["cds"] = sorted(segs)
    for tid, segs in exons.items():
        if tid in txs:
            txs[tid]["exons"] = sorted(segs)
    log("注释解析：%d 行 / %d 转录本 / CDS %d 段 / exon %d 段"
        % (stats["lines"], len(txs), stats["cds"], stats["exon"]))
    return txs, stats


# ---------------------------------------------------------------- 建层
def splice(tx, seqs, log=None):
    """按转录本坐标拼接 CDS，返回 (nt, 连接点 nt 偏移列表, 问题列表)。"""
    problems = []
    segs = tx["cds"]
    if not segs:
        return "", [], ["no_cds"]
    chrom = seqs.get(tx["seqid"])
    if chrom is None:
        return "", [], ["seqid_missing:%s" % tx["seqid"]]
    order = sorted(segs, key=lambda s: s[0], reverse=(tx["strand"] == "-"))
    parts, cuts, acc = [], [], 0
    for idx, (s, e, _ph) in enumerate(order):
        if s < 1 or e > len(chrom) or e < s:
            problems.append("bad_coord:%d-%d" % (s, e))
            return "", [], problems
        sub = chrom[s - 1:e]
        if tx["strand"] == "-":
            sub = revcomp(sub)
        parts.append(sub)
        if idx < len(order) - 1:
            acc += len(sub)
            cuts.append(acc)
    return "".join(parts), cuts, problems


def junctions(nt, cuts, flank_aa=12, table=1):
    """从拼接后的 nt 造跨连接点肽（两侧各 flank_aa 个 aa）。

    只接受密码子对齐的连接点（cut % 3 == 0）且两侧无终止密码子。
    """
    out = []
    for cut in cuts:
        if cut % 3:
            continue
        left = nt[max(0, cut - 3 * flank_aa):cut]
        right = nt[cut:cut + 3 * flank_aa]
        if len(left) < 3 or len(right) < 3:
            continue
        pep = translate(left + right, table)
        if "*" in pep or len(pep) < 5:
            continue
        out.append((cut, pep))
    return out


def tx_orfs(seq, min_aa=30, table=1):
    """转录本 3 框 ORF（按终止密码子切段），返回 [(frame, start_nt, pep)]。"""
    out = []
    codons = CODONS.get(table) or CODONS[1]
    for frame in range(3):
        nt = seq[frame:]
        start = 0
        for i in range(0, len(nt) - len(nt) % 3, 3):
            if codons.get(nt[i:i + 3].upper()) == "*":
                pep = translate(nt[start:i], table)
                if len(pep) >= min_aa:
                    out.append((frame + 1, start, pep))
                start = i + 3
        pep = translate(nt[start:], table, stop_at_stop=True)
        if len(pep) >= min_aa:
            out.append((frame + 1, start, pep))
    return out


def build(workdir, sp, fna, gff, transcripts=None, flank_aa=12, min_aa=30,
          ref_protein=None, preset=None, log=print):
    """建真核层；返回统计字典（同时写 manifest 与 stats.tsv）。"""
    out_dir = os.path.join(workdir, "02_db", sp)
    os.makedirs(out_dir, exist_ok=True)
    seqs = {k: v.upper() for k, v in read_fasta(fna).items()}
    log("基因组：%d 条序列 / %.1f Mb"
        % (len(seqs), sum(len(v) for v in seqs.values()) / 1e6))
    txs, parse_stats = parse_annotation(gff, log=log)

    prot_path = os.path.join(out_dir, "L0_annot.faa")
    junc_path = os.path.join(out_dir, "L0_junction.faa")
    rows = []
    n_prot = 0
    junc_seen = {}
    junc_rows = []
    reasons = defaultdict(int)
    with open(prot_path, "w", encoding="utf-8", newline="\n") as fh:
        for tid in sorted(txs):
            tx = txs[tid]
            nt, cuts, problems = splice(tx, seqs)
            if problems:
                reasons[problems[0].split(":")[0]] += 1
                rows.append([tid, tx.get("gene", ""), tx["seqid"], tx["strand"],
                             0, 0, 0, 0, problems[0]])
                continue
            table = tx.get("table", 1)
            pep = translate(nt, table, stop_at_stop=True)
            # 真正的内部终止：终止符出现在最后一个密码子之前（末尾的终止不算）
            if "*" in translate(nt, table)[:-1]:
                reasons["internal_stop"] += 1
            jp = junctions(nt, cuts, flank_aa, table)
            for k, (cut, p) in enumerate(jp):
                junc_rows.append((tid, cut, p))
                junc_seen.setdefault(p, 0)
                junc_seen[p] += 1
            if len(pep) < min_aa:
                reasons["short"] += 1
                rows.append([tid, tx.get("gene", ""), tx["seqid"], tx["strand"],
                             len(tx["cds"]), len(nt), len(pep), len(jp), "short"])
                continue
            n_prot += 1
            fh.write(">L0|%s|gene=%s|chr=%s|strand=%s|exons=%d|len=%d\n%s\n"
                     % (tid, tx.get("gene", ""), tx["seqid"], tx["strand"],
                        len(tx["cds"]), len(pep), pep))
            rows.append([tid, tx.get("gene", ""), tx["seqid"], tx["strand"],
                         len(tx["cds"]), len(nt), len(pep), len(jp), "ok"])

    with open(junc_path, "w", encoding="utf-8", newline="\n") as fh:
        n_junc = 0
        for tid, cut, pep in junc_rows:
            n_junc += 1
            fh.write(">LJR|%s|nt=%d|len=%d\n%s\n" % (tid, cut, len(pep), pep))

    n_txorf = 0
    txorf_path = os.path.join(out_dir, "L2_txorf.faa")
    if transcripts:
        txseqs = read_fasta(transcripts)
        with open(txorf_path, "w", encoding="utf-8", newline="\n") as fh:
            for tid in sorted(txseqs):
                for frame, start, pep in tx_orfs(txseqs[tid], min_aa=min_aa):
                    n_txorf += 1
                    fh.write(">L2|%s|f%d|nt=%d|len=%d\n%s\n"
                             % (tid, frame, start, len(pep), pep))
        log("转录本层：%d 条转录本 → %d 个 ORF" % (len(txseqs), n_txorf))

    stats_path = os.path.join(out_dir, "euk_stats.tsv")
    with open(stats_path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("transcript\tgene\tchr\tstrand\tn_cds\tnt_len\taa_len\t"
                 "junctions\tstatus\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")

    manifest = {
        "species": sp,
        "inputs": {"fna": os.path.abspath(fna), "gff": os.path.abspath(gff),
                   "transcripts": os.path.abspath(transcripts) if transcripts else None},
        "params": {"flank_aa": flank_aa, "min_aa": min_aa},
        "counts": {
            "genome_seqs": len(seqs),
            "genome_mb": round(sum(len(v) for v in seqs.values()) / 1e6, 2),
            "transcripts": len(txs),
            "proteins": n_prot,
            "junction_peptides": n_junc,
            "junction_unique": len(junc_seen),
            "tx_orfs": n_txorf,
            "skipped": dict(reasons),
        },
        "outputs": {"L0_annot.faa": n_prot,
                    "L0_junction.faa": n_junc,
                    "L2_txorf.faa": n_txorf},
    }
    if ref_protein:
        manifest["reference_check"] = compare_with_reference(prot_path, ref_protein)
    if preset:
        from . import eukpresets
        chk = eukpresets.check(n_prot, preset)
        manifest["preset"] = {"code": preset, "check": chk}
        if chk.get("expected"):
            log("预设核对（%s）：蛋白 %d / 文档 %d，偏差 %+.2f%%"
                % (preset, n_prot, chk["expected"], chk["delta_pct"]))
    with open(os.path.join(out_dir, "euk_manifest.json"), "w",
              encoding="utf-8") as fh:
        json.dump(manifest, fh, ensure_ascii=False, indent=1)
    log("完成：蛋白 %d / 连接肽 %d（唯一 %d）/ 转录本 ORF %d → %s"
        % (n_prot, n_junc, len(junc_seen), n_txorf, out_dir))
    return manifest


_IL = str.maketrans("IL", "LL")


def il_norm(seq):
    """I→L 归一（pFind 判肽口径；L0 去重也用它）。"""
    return (seq or "").upper().translate(_IL)


def il_dedupe_stats(path):
    """统计 I/L 归一去重后的蛋白数：{records, unique}。

    人睾丸文档的口径：RefSeq 136,807 → 90,373 条唯一（+CON 571 = db0 90,944）。
    """
    recs = read_fasta(path)
    seen = set()
    for seq in recs.values():
        seen.add(il_norm(seq))
    return {"records": len(recs), "unique": len(seen)}


def compare_with_reference(our_faa, ref_faa, near_pct=0.02):
    """与参考蛋白集比对：exact / prefix（末尾少 1 个密码子）/ near（少数替换）。

    只比序列集合（参考头格式各异，按 ID 对齐容易错）。`near` 用前 40 aa 做键，
    命中的再算 Hamming——人类这类"RefSeq 与基因组有若干碱基差异"的情况会落在
    这里（GFF 的 Note 里就写着 "has N substitutions compared to this genomic
    sequence"），不算剪接错。
    """
    ours = read_fasta(our_faa)
    refs = set(read_fasta(ref_faa).values())
    by_len = defaultdict(list)
    by_prefix = defaultdict(list)
    for r in refs:
        by_len[len(r)].append(r)
        if len(r) >= 40:
            by_prefix[r[:40]].append(r)
    exact = prefix = near = 0
    for pep in ours.values():
        if pep in refs:
            exact += 1
            continue
        short = []
        for L in range(len(pep) - 3, len(pep) + 4):
            short.extend(by_len.get(L, ()))
        if any(r.startswith(pep) or pep.startswith(r) for r in short):
            prefix += 1
            continue
        cands = [r for r in by_prefix.get(pep[:40], ()) if len(r) == len(pep)]
        if cands:
            best = min(sum(1 for a, b in zip(pep, r) if a != b) for r in cands)
            if best <= max(2, int(len(pep) * near_pct)):
                near += 1
    total = len(ours) or 1
    return {"ours": len(ours), "reference": len(refs), "exact": exact,
            "prefix_1aa": prefix, "near_sub": near,
            "exact_rate": round(exact / total, 4),
            "match_rate": round((exact + prefix + near) / total, 4)}
