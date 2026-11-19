# -*- coding: utf-8 -*-
"""基因组六框翻译是否适用：按基因组紧凑度给判定，而不是拍脑袋。

判据只有两个数：
  d = ΣCDS / 基因组  —— 编码密度
  R = 六框库残基 / 注释蛋白残基 ≈ 2G / P ≈ 6 / d
     （六框共 6 个读框，每框产出 1/3 个残基/nt → 总残基 ≈ 2G；
       注释蛋白 P ≈ ΣCDS/3，所以 R 与 d 一一对应）

R 的物理含义：六框库比"已知蛋白"大多少倍。库越大，同一 FDR 下的分数门槛越高、
真肽越容易被挤掉，所以 R 就是六框该不该做的直接度量。

分级（阈值可按物种群调整，写在 TIERS 里）：
  compact   R ≤ 15（d ≥ 40%）：推荐六框，库与转录本层同一量级
  moderate  15 < R ≤ 30（20% ≤ d < 40%）：无转录组时才做
  sparse    30 < R ≤ 60（10% ≤ d < 20%）：不建议
  extreme   R > 60（d < 10%）：不做（改转录本/组装 + 从头预测）
另有绝对算力闸：六框残基 > 150 M aa 或 ORF > 2 M 时，即便紧凑也不做。

用法：
    python -m orfscape sframe --workdir <proj> --sp H99
    python -m orfscape sframe --workdir <proj> --sp H99 --no-measure   # 只算解析量，不真翻译
"""
import json
import os

# 判定阈值（R = 六框残基 / 注释蛋白残基）
TIERS = (
    ("compact", 15.0, "推荐做六框", "six-frame recommended"),
    ("moderate", 30.0, "仅在无转录组时做六框", "six-frame only if no transcriptome"),
    ("sparse", 60.0, "不建议做六框", "six-frame not advised"),
    ("extreme", 1e9, "不做六框", "do not build a six-frame DB"),
)

ABS_RES_MAX = 150_000_000      # 六框残基算力闸
ABS_ORF_MAX = 2_000_000        # ORF 条数算力闸

_BASES = "TCAG"
_AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
_TABLE = {}
_i = 0
for _a in _BASES:
    for _b in _BASES:
        for _c in _BASES:
            _TABLE[_a + _b + _c] = _AAS[_i]
            _i += 1
_COMP = str.maketrans("ACGTacgt", "TGCAtgca")


def read_fasta_len(path):
    """(条数, 总 nt)——只数字符，不建序列对象，大基因组也不爆内存。"""
    n = nt = 0
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                n += 1
            else:
                nt += len(line.strip())
    return n, nt


def genome_metrics(root, sp):
    """基因组 / 注释 / 蛋白 / 转录本四个量。"""
    gdir = os.path.join(str(root), "01_genome", sp)
    out = {"species": sp}
    fa = os.path.join(gdir, "genomic.fna")
    if not os.path.isfile(fa):
        raise FileNotFoundError(fa)
    out["contigs"], out["genome_bp"] = read_fasta_len(fa)

    cds_bp = cds_n = genes = 0
    gff = os.path.join(gdir, "genomic.gff")
    if os.path.isfile(gff):
        with open(gff, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 9:
                    continue
                try:
                    if f[2] == "CDS":
                        cds_n += 1
                        cds_bp += int(f[4]) - int(f[3]) + 1
                    elif f[2] == "gene":
                        genes += 1
                except ValueError:
                    continue
    out["cds_n"], out["cds_bp"], out["genes"] = cds_n, cds_bp, genes

    prot = os.path.join(gdir, "protein.faa")
    pn = pres = 0
    if os.path.isfile(prot):
        pn, pres = read_fasta_len(prot)
    out["proteins"], out["protein_res"] = pn, pres

    tx_nt = tx_n = 0
    tdir = os.path.join(str(root), "02_tx", sp)
    for name in ("transcripts.fa", "transcripts_aug.fa"):
        p = os.path.join(tdir, name)
        if os.path.isfile(p):
            n, nt = read_fasta_len(p)
            tx_n += n
            tx_nt += nt
            out[name] = {"n": n, "nt": nt}
    out["tx_n"], out["tx_nt"] = tx_n, tx_nt
    return out


def _translate(seq):
    s = seq.upper().replace("N", "A")
    return "".join(_TABLE.get(s[i:i + 3], "X") for i in range(0, len(s) - 2, 3))


def _load_gene_spans(gff):
    """参考注释的基因区间（1-based 闭区间），按 contig 排序。"""
    import bisect  # noqa: F401  （保留给调用方按需二分）
    spans = {}
    with open(gff, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            try:
                spans.setdefault(f[0], []).append((int(f[3]), int(f[4])))
            except ValueError:
                continue
    for k in spans:
        spans[k].sort()
    return spans


def _overlaps(spans, contig, lo, hi):
    a = spans.get(contig, [])
    if not a:
        return False
    import bisect
    i = bisect.bisect_right([x[0] for x in a], hi)
    for j in range(max(0, i - 3), min(len(a), i + 3)):
        if a[j][0] <= hi and a[j][1] >= lo:
            return True
    return False


def sixframe_measure(genome_fna, min_aa=30, gene_gff=None):
    """真翻译一遍：六框 ORF 条数、残基、含 X 条数，以及基因内/基因间拆分。

    基因间占比决定了"整库六框"值不值得：占绝大多数 ORF 都在已注释基因内的
    不同读框上时，整库只是把库放大，真正新增的假设很有限。
    """
    spans = _load_gene_spans(gene_gff) if gene_gff else {}
    seqs, cur = {}, None
    with open(genome_fna, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                cur = line[1:].split()[0]
                seqs[cur] = []
            elif cur:
                seqs[cur].append(line)
    n = res = with_x = ig = ig_res = 0
    for name, parts in seqs.items():
        s = "".join(parts)
        L = len(s)
        for strand, seq in (("+", s), ("-", s.translate(_COMP)[::-1])):
            for frame in range(3):
                prot = _translate(seq[frame:])
                off = frame
                for orf in prot.split("*"):
                    ln = len(orf)
                    lo, hi = off, off + ln * 3 - 1
                    off += ln * 3 + 3
                    if ln < min_aa:
                        continue
                    n += 1
                    if "X" in orf:
                        with_x += 1
                        continue
                    res += ln
                    if spans:
                        g0, g1 = ((lo, hi) if strand == "+"
                                  else (L - 1 - hi, L - 1 - lo))
                        if not _overlaps(spans, name, g0, g1):
                            ig += 1
                            ig_res += ln
    return {"sixframe_orfs": n, "sixframe_res": res, "orfs_with_x": with_x,
            "sixframe_intergenic_orfs": ig, "sixframe_intergenic_res": ig_res}


def verdict(m, min_aa=30):
    """按 R（= 六框/注释蛋白）与绝对算力闸给结论。"""
    G, P = m.get("genome_bp") or 0, m.get("protein_res") or 0
    d = (m.get("cds_bp") or 0) / float(G) if G else 0.0
    R = (2.0 * G / P) if P else None
    if not R:
        level, zh, en = "unknown", "缺注释蛋白，无法判定", "no annotated proteome"
    else:
        level, zh, en = "extreme", "不做六框", "do not build a six-frame DB"
        for name, thr, z, e in TIERS:
            if R <= thr:
                level, zh, en = name, z, e
                break
    res6 = m.get("sixframe_res") or int(2 * G)
    orfs6 = m.get("sixframe_orfs")
    block = []
    if res6 > ABS_RES_MAX:
        block.append("六框残基 %s > %s（算力）"
                     % ("{:,}".format(res6), "{:,}".format(ABS_RES_MAX)))
    if orfs6 and orfs6 > ABS_ORF_MAX:
        block.append("六框 ORF %s > %s（算力）"
                     % ("{:,}".format(orfs6), "{:,}".format(ABS_ORF_MAX)))
    if block:
        level, zh, en = "extreme", "不做六框（" + "；".join(block) + "）", "do not build"
    out = {"level": level, "verdict_zh": zh, "verdict_en": en,
           "coding_density": round(d, 4), "inflation_R": (round(R, 2) if R else None),
           "blocked_by": block}
    # 基因间占比：决定"整库六框"还是"只做基因间六框"
    ig = m.get("sixframe_intergenic_orfs")
    tot = m.get("sixframe_orfs")
    if ig is not None and tot:
        frac = ig / float(tot)
        out["intergenic_fraction"] = round(frac, 4)
        if level == "compact" and frac < 0.15:
            out["scope_zh"] = ("建议只建**基因间**六框：整库 %s 条里只有 %s 条（%.1f%%）"
                              "落在注释基因之外，其余都是已注释基因内的不同读框"
                              % ("{:,}".format(tot), "{:,}".format(ig), 100 * frac))
            out["scope_en"] = ("build intergenic-only six-frame: only %s of %s ORFs "
                               "(%.1f%%) lie outside annotated genes"
                               % ("{:,}".format(ig), "{:,}".format(tot), 100 * frac))
        elif level == "compact":
            out["scope_zh"] = "可整库六框（基因间占 %.1f%%）" % (100 * frac)
            out["scope_en"] = "whole-genome six-frame acceptable (%.1f%% intergenic)" % (
                100 * frac)
    return out


def run(root, sp, min_aa=30, measure=True, log=print):
    m = genome_metrics(root, sp)
    if measure:
        gdir = os.path.join(str(root), "01_genome", sp)
        m.update(sixframe_measure(os.path.join(gdir, "genomic.fna"),
                                  min_aa=min_aa,
                                  gene_gff=os.path.join(gdir, "genomic.gff")))
    v = verdict(m, min_aa=min_aa)
    m.update(v)

    log("  [%s] 基因组 %s bp（%d contig）· CDS %s bp / %s 条 · 编码密度 %.1f%%"
        % (sp, "{:,}".format(m["genome_bp"]), m["contigs"],
           "{:,}".format(m["cds_bp"]), "{:,}".format(m["cds_n"]),
           100 * v["coding_density"]))
    log("  注释蛋白 %s 条 / %s 残基"
        % ("{:,}".format(m["proteins"]), "{:,}".format(m["protein_res"])))
    if measure:
        log("  六框实测（≥%d aa）：ORF %s 条 / 残基 %s（含 X 的 %s 条已剔除）"
            % (min_aa, "{:,}".format(m["sixframe_orfs"]),
               "{:,}".format(m["sixframe_res"]), "{:,}".format(m["orfs_with_x"])))
        if m.get("sixframe_orfs"):
            log("  其中落在注释基因之外（基因间）：%s 条（%.1f%%）/ 残基 %s"
                % ("{:,}".format(m.get("sixframe_intergenic_orfs", 0)),
                   100.0 * m.get("sixframe_intergenic_orfs", 0) / m["sixframe_orfs"],
                   "{:,}".format(m.get("sixframe_intergenic_res", 0))))
    else:
        log("  六框估算：残基 ≈ 2G = %s（未真翻译）" % "{:,}".format(2 * m["genome_bp"]))
    if m.get("inflation_R"):
        log("  膨胀倍数 R = 六框/注释蛋白 = %.1f×（= 6/编码密度 %.1f×）"
            % (m["inflation_R"], 6.0 / max(v["coding_density"], 1e-9)))
    if m.get("tx_nt"):
        log("  对照：转录本 %s 条 / %s nt → 三框残基 ≈ %s（%.1f M）"
            % ("{:,}".format(m["tx_n"]), "{:,}".format(m["tx_nt"]),
               "{:,}".format(m["tx_nt"]), m["tx_nt"] / 1e6))
        if measure and m["tx_nt"]:
            log("  六框 / 转录本层 = %.2f×（≤1.5 时六框几乎不增加库规模）"
                % (m["sixframe_res"] / float(m["tx_nt"])))
    log("  判定：%s" % v["verdict_zh"])
    if m.get("scope_zh"):
        log("  范围：%s" % m["scope_zh"])
    return m


def write_plan(d, path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(d, fh, ensure_ascii=False, indent=2, sort_keys=True)
        fh.write("\n")
    return path
