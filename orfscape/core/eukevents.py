# -*- coding: utf-8 -*-
"""真核新肽事件分类（基因间 / N 端 / C 端 / 内含子 / 不同框 / 同框 / 反义 / ncRNA）与肽级归一。

口径参考 pAnno v3（2026-11-11 裁定）并适配到 ORFscape 的层命名：

  1. ncRNA / lncRNA 基因内的肽 → ``ncRNA_translation``（不再并进 intronic）；
  2. 延伸判定**以肽段坐标为准**：肽段（或 ORF）落在注释 CDS 外显子之外才算
     ``nterm_ext`` / ``cterm_ext``；ORF 几何上超出 CDS、但肽段全在 CDS 内 →
     ``alt_frame``；
  3. 跨外显子连接点的肽（``LJR|`` 层）按**分段基因组区间**判定，不按连续区间；
  4. 事件计数**肽级归一**：每条肽只保留一个主事件（Tier 最优优先，其次有基因归属），
     避免 ORF 端计数把同一肽重复计入。

支持的蛋白名（``eukbuild`` 产出）：

  ``L0|<tid>|gene=..|chr=..|strand=..|exons=..|len=..``   注释（剪接）蛋白
  ``LJR|<tid>|nt=<cut>|len=<n>``                          剪接连接肽
  ``L2|<tid>|f<frame>|nt=<start>|len=<n>``                转录本 3 框 ORF

坐标系：GFF 1-based 闭区间；内部一律换成 0-based 闭区间再比较。
"""
import collections
import os
import re
from pathlib import Path

from .novel import iln, load_spectra
from .report import load_tsv

# 七类 + 一个"注释内同框"兜底
INTERGENIC = "intergenic"
NTERM = "nterm_ext"
CTERM = "cterm_ext"
INTRONIC = "intronic"
ALT_FRAME = "alt_frame"
ANTISENSE = "antisense"
NCRNA = "ncRNA_translation"
INTERNAL = "inframe_internal"

CLASSES = (INTERGENIC, NTERM, CTERM, INTRONIC, ALT_FRAME, ANTISENSE, NCRNA,
           INTERNAL)

CLASS_LABEL = {
    INTERGENIC: ("基因间区", "Intergenic"),
    NTERM: ("N 端延伸", "N-terminal extension"),
    CTERM: ("C 端延伸", "C-terminal extension"),
    INTRONIC: ("内含子内", "Intronic"),
    ALT_FRAME: ("不同读框", "Alternative frame"),
    ANTISENSE: ("反义链", "Antisense"),
    NCRNA: ("非编码 RNA 翻译", "ncRNA translation"),
    INTERNAL: ("注释 CDS 内同框", "In-frame inside annotated CDS"),
}

# Ig / TCR 基因座（人、免疫组织必须排除；pAnno v3 同口径）
_IG_TCR = re.compile(r"^(IGH|IGK|IGL|TCR)[A-Z0-9]*$", re.I)


def _open_w(base, path, mode="w", **kw):
    """写盘前确认目标落在 base 之下：输出路径都来自 CLI，必须先收敛。"""
    p = os.path.abspath(str(path))
    b = os.path.abspath(str(base))
    if os.path.commonpath([p, b]) != b:
        raise ValueError("路径越界: %s 不在 %s 之下" % (p, b))
    return open(p, mode, **kw)


def _is_ig_tcr(name):
    return bool(_IG_TCR.match(str(name or "").split(".")[0]))


def _seg(name, fallback=""):
    """单段路径名（物种代码 / 子目录）：拒绝分隔符与父目录跳转。"""
    s = str(name or "").strip()
    if (not s or s in (".", "..") or "/" in s or "\\" in s or os.path.isabs(s)):
        return fallback
    return s


def _within(base, *parts):
    """base 之下的路径：规范化后必须仍在 base 内，否则抛 ValueError。"""
    root0 = os.path.realpath(base)
    p = os.path.realpath(os.path.join(root0, *parts))
    if p != root0 and not p.startswith(root0 + os.sep):
        raise ValueError("path escapes base: %s" % p)
    return p


def parse_euk_name(name):
    """``L0|`` / ``LJR|`` / ``L2|`` 蛋白名 → (kind, tid, extra dict) 或 None。"""
    f = str(name or "").split("|")
    if len(f) < 2 or f[0] not in ("L0", "LJR", "L2"):
        return None
    extra = {}
    for token in f[2:]:
        if "=" in token:
            k, v = token.split("=", 1)
            extra[k] = v
        elif len(token) > 1 and token[0] == "f" and token[1:].isdigit():
            extra["f"] = token          # L2|<tid>|f<frame> 的读框位
    return f[0], f[1], extra


# ---------------------------------------------------------------- 转录本模型
def load_tx_model(gff, log=None):
    """GFF → (txs, genes)。

    txs   : {tid: {seqid, strand, gene, biotype, exons[(s,e)], cds[(s,e,phase)], nt_len}}
    genes : {gid: {seqid, strand, biotype, exons, cds, lo, hi}}    (1-based 闭区间)
    """
    from .eukbuild import parse_annotation

    gff = os.path.realpath(str(gff))
    if not os.path.isfile(gff):
        raise FileNotFoundError("GFF not found: %s" % gff)
    txs, _stats = parse_annotation(gff, log=lambda *_: None)
    genes = {}
    for tid, tx in txs.items():
        gid = tx.get("gene") or tid
        g = genes.setdefault(gid, {"seqid": tx["seqid"], "strand": tx["strand"],
                                   "biotype": tx.get("biotype", ""),
                                   "exons": [], "cds": []})
        g["exons"].extend(tx.get("exons") or [])
        g["cds"].extend(tx.get("cds") or [])
    out_txs = {}
    for tid, tx in txs.items():
        exons = tx.get("exons") or [(s, e) for s, e, _p in (tx.get("cds") or [])]
        exons = sorted(exons)
        nt_len = sum(e - s + 1 for s, e in exons)
        out_txs[tid] = {"seqid": tx["seqid"], "strand": tx["strand"],
                        "gene": tx.get("gene", ""), "biotype": tx.get("biotype", ""),
                        "exons": exons, "cds": list(tx.get("cds") or []),
                        "nt_len": nt_len}
    for g in genes.values():
        if not g["exons"]:
            g["exons"] = [(s, e) for s, e, _p in g["cds"]]
        g["exons"].sort()
        g["lo"] = min((s for s, _e in g["exons"]), default=0)
        g["hi"] = max((e for _s, e in g["exons"]), default=0)
        g["cds"].sort()
    if log:
        log("  真核注释：%d 转录本 / %d 基因" % (len(out_txs), len(genes)))
    return out_txs, genes


def tx_blocks(tx, nt_start, nt_len):
    """转录本 nt 区间 → 分段基因组区间（0-based 闭区间，按链方向走外显子）。"""
    exons = tx.get("exons") or []
    if not exons or nt_len <= 0:
        return []
    order = sorted(exons) if tx["strand"] == "+" else sorted(exons, reverse=True)
    blocks, pos = [], 0
    end = nt_start + nt_len
    for s, e in order:
        length = e - s + 1
        a, b = max(nt_start, pos), min(end, pos + length)
        if a < b:
            if tx["strand"] == "+":
                lo1, hi1 = s + (a - pos), s + (b - pos) - 1
            else:
                lo1, hi1 = e - (b - pos - 1), e - (a - pos)
            blocks.append((tx["seqid"], tx["strand"], lo1 - 1, hi1 - 1))
        pos += length
    return blocks


_L1SF_RX = re.compile(r"^L1_SF\|([^|]+)\|([+-])f(\d)\|nt=(\d+)\|len=(\d+)$")


def _l1sf(name):
    """`L1_SF|<contig>|<strand>f<frame>|nt=<0基起点>|len=<aa>` → 四元组。

    基因间区六框的 ORF 是基因组上连续一段，不需要转录本模型：
    起点 nt 是该 ORF 第一个碱基在**正向基因组**上的 0 基坐标（负链为其 5′ 端）。
    """
    m = _L1SF_RX.match(str(name or ""))
    if not m:
        return None
    return m.group(1), m.group(2), int(m.group(4)), int(m.group(5))


def pep_blocks(name, pos1, pep_len, txs):
    """蛋白名 + pFind 位点 → 肽段的分段基因组区间（0-based 闭区间）。

    pos1 是 pFind 的 0-based 氨基酸起点（见口径：aa→nt 一律 pos*3）。
    """
    sf = _l1sf(name)
    if sf:                                  # L1_SF：连续 ORF，直接算
        contig, strand, nt, _aa = sf
        pos0 = int(pos1)
        if strand == "+":
            lo = nt + pos0 * 3
            return [(contig, "+", lo, lo + pep_len * 3 - 1)]
        hi = nt - pos0 * 3
        return [(contig, "-", hi - pep_len * 3 + 1, hi)]
    parsed = parse_euk_name(name)
    if not parsed:
        return []
    kind, tid, extra = parsed
    tx = txs.get(tid)
    if not tx:
        return []
    pos0 = int(pos1)
    if kind == "L0":
        start_nt = pos0 * 3
    elif kind == "LJR":
        try:
            cut = int(extra.get("nt", "0"))
        except ValueError:
            return []
        right = min(36, max(tx["nt_len"] - cut, 0))
        left = max(3 * pep_len - right, 0)
        start_nt = max(cut - left, 0) + pos0 * 3
    else:                                        # L2|f<frame>|nt=<start>
        try:
            frame = int(str(extra.get("f", "1")).lstrip("f"))
            start = int(extra.get("nt", "0"))
        except ValueError:
            return []
        start_nt = (frame - 1) + start + pos0 * 3
    return tx_blocks(tx, start_nt, pep_len * 3)


def orf_blocks(name, txs):
    """ORF 自身的基因组跨度（用于判"ORF 延伸但肽段没延伸"）。"""
    sf = _l1sf(name)
    if sf:
        contig, strand, nt, aa = sf
        if strand == "+":
            return [(contig, "+", nt, nt + aa * 3 - 1)]
        return [(contig, "-", nt - aa * 3 + 1, nt)]
    parsed = parse_euk_name(name)
    if not parsed:
        return []
    kind, tid, extra = parsed
    tx = txs.get(tid)
    if not tx:
        return []
    try:
        length = int(extra.get("len", "0"))
    except ValueError:
        return []
    if kind == "LJR":
        # 连接肽没有单独的 start 字段：按"以连接点为中心、总长 len"近似
        try:
            cut = int(extra.get("nt", "0") or 0)
        except ValueError:
            return []
        return tx_blocks(tx, max(cut - 3 * (length // 2), 0), length * 3)
    if kind == "L0":
        # 注释蛋白本身就是 CDS 的翻译产物，没有"ORF 延伸"可言
        return []
    try:
        frame = int(str(extra.get("f", "1")).lstrip("f"))
        start = int(extra.get("nt", "0"))
    except ValueError:
        return []
    return tx_blocks(tx, (frame - 1) + start, length * 3)


# ---------------------------------------------------------------- 分类
def _overlap(lo1, hi1, lo2, hi2):
    return not (hi1 < lo2 or hi2 < lo1)


def _is_ncrna(gene):
    return (not gene["cds"]) or ("RNA" in (gene.get("biotype") or "").upper()
                                 and "MRNA" not in (gene.get("biotype") or "").upper())


def _gene_of_overlap(blocks, genes_by_contig, strand):
    """返回 (same_strand_hits, antisense_hits)，按重叠长度降序。"""
    same, anti = [], []
    seen = set()
    for contig, _st, lo, hi in blocks:
        for g in genes_by_contig.get(contig, []):
            key = (g["seqid"], g["lo"], g["hi"], g["strand"], g.get("gene", ""))
            if key in seen:
                continue
            if not _overlap(lo, hi, g["lo"] - 1, g["hi"] - 1):
                continue
            seen.add(key)
            (same if g.get("strand") == strand else anti).append(g)
    same.sort(key=lambda g: (g["hi"] - g["lo"]))
    return same, anti


def _in_frame(pep_blocks_, cds, strand):
    """肽段分段区间与注释 CDS 是否同框（按重叠 CDS 的相位算）。"""
    for _contig, _st, lo, hi in pep_blocks_:
        for s, e, ph in cds:
            if not _overlap(lo, hi, s - 1, e - 1):
                continue
            if strand == "+":
                rel = lo - (s - 1) - ph
            else:
                rel = (e - 1) - ph - hi
            if rel % 3 == 0:
                return True
    return False


def _min_gap(blocks, genes_by_contig):
    """肽段到最近注释基因跨度的距离（bp）；0 = 重叠，None = 该 contig 无基因。

    基因间/基因内边界处的肽用这个数区分"基因沙漠里的新 ORF"与"贴着基因边界的延伸"。
    """
    best = None
    for contig, _st, lo, hi in blocks:
        for g in genes_by_contig.get(contig, []):
            glo, ghi = g["lo"] - 1, g["hi"] - 1
            if _overlap(lo, hi, glo, ghi):
                d = 0
            elif ghi < lo:          # 基因在左：相邻也算 1 bp，避免与"重叠"混成 0
                d = lo - ghi
            else:                   # 基因在右
                d = glo - hi
            if best is None or d < best:
                best = d
    return best


def classify(pep_blocks_, orf_spans, gene, strand):
    """单条肽在单基因内的分类（肽段坐标驱动）。"""
    cds = gene["cds"]
    cds_lo = min(s for s, _e, _p in cds)
    cds_hi = max(e for _s, e, _p in cds)
    cds0 = [(s - 1, e - 1) for s, e, _p in cds]

    def inside_any_cds(block):
        return any(_overlap(block[2], block[3], a, b) for a, b in cds0)

    outside = [b for b in pep_blocks_ if not inside_any_cds(b)]
    if outside:
        n_side = c_side = False
        for _contig, _st, lo, hi in outside:
            if strand == "+":
                n_side = n_side or lo < cds_lo - 1
                c_side = c_side or hi > cds_hi - 1
            else:
                n_side = n_side or hi > cds_hi - 1
                c_side = c_side or lo < cds_lo - 1
        if n_side and not c_side:
            return NTERM
        if c_side and not n_side:
            return CTERM
        if n_side and c_side:
            return NTERM                   # 双侧都延伸：按 pAnno 取 N 侧
        return INTRONIC                    # 肽段在 CDS 之外、但夹在首末 CDS 之间（内含子区）
    # 肽段全在 CDS 内：ORF 是否几何延伸？
    orf_outside = False
    for _contig, _st, lo, hi in orf_spans:
        if strand == "+":
            orf_outside = orf_outside or lo < cds_lo - 1 or hi > cds_hi - 1
        else:
            orf_outside = orf_outside or hi > cds_hi - 1 or lo < cds_lo - 1
    if orf_outside:
        return ALT_FRAME                   # ORF 延伸但肽段无延伸证据
    return INTERNAL if _in_frame(pep_blocks_, cds, strand) else ALT_FRAME


def classify_peptide(name, pos1, pep_len, txs, genes_by_contig):
    """一条肽 → (class, gene_symbol)。"""
    blocks = pep_blocks(name, pos1, pep_len, txs)
    if not blocks:
        return INTERGENIC, ""
    strand = blocks[0][1]
    same, anti = _gene_of_overlap(blocks, genes_by_contig, strand)
    same = [g for g in same if not _is_ig_tcr(g.get("gene"))]
    if same:
        gene = same[0]
        if _is_ncrna(gene):
            return NCRNA, gene.get("gene", "")
        return classify(blocks, orf_blocks(name, txs), gene, strand), gene.get("gene", "")
    if anti:
        return ANTISENSE, anti[0].get("gene", "")
    return INTERGENIC, ""


# ---------------------------------------------------------------- 肽级归一
_TIER_RANK = {"T1": 0, "T2": 1, "T3": 2, "Tier1": 0, "Tier2": 1, "Tier3": 2}


def normalize_peptide_level(hits):
    """{peptide: [hit, ...]} → 每肽一个主事件。

    hit 形如 ``{class, gene, tier, orf, pos}``；排序键 = (Tier 优先, 有基因归属优先)。
    """
    main = {}
    for pep, items in hits.items():
        best = sorted(items, key=lambda h: (_TIER_RANK.get(h.get("tier"), 2),
                                            0 if h.get("gene") else 1))[0]
        classes = sorted({h["class"] for h in items})
        main[pep] = dict(best)
        main[pep]["all_classes"] = classes
        main[pep]["n_events"] = len(items)
    return main


def summarize(main):
    counts = collections.Counter(v["class"] for v in main.values())
    return {k: counts.get(k, 0) for k in CLASSES if counts.get(k, 0)}


# ---------------------------------------------------------------- 入口
def run(workdir, sp, rows, out_dir=None, log=print):
    """P7 可信肽 → 真核事件分类表。

    rows 需要 ``peptide`` / ``tier`` / 以及 ``protein``+``pos0`` 或坐标表里的位点；
    返回 {peptide: 主事件 dict}。
    """
    # 转录本模型优先级：02_tx/<sp>/tx_model.gff（原生注释 + 外部组装新增转录本）
    # > 01_genome/<sp>/genomic.gff。L2| 的 tid 可能是外源组装的，只有前者能映射。
    gff = os.path.join(workdir, "02_tx", _seg(sp), "tx_model.gff")
    if not os.path.isfile(gff):
        gff = os.path.join(workdir, "01_genome", _seg(sp), "genomic.gff")
    if not os.path.isfile(gff):
        log("  缺 %s，真核事件分类跳过" % gff)
        return {}
    # 坐标映射用 tx_model（含外部组装/AUGUSTUS 转录本，L2| 的 tid 只在里面）；
    # **基因几何只用参考注释**：tx_model 里 MSTRG/AUG 转录本按区间挂到注释基因上，
    # 会把基因 span 撑到基因间区，让 intergenic 恒为 0；而没有 CDS 的新模型
    # 又会被当成"非编码基因"，把基因间的新 ORF 误判成 ncRNA_translation。
    txs, _model_genes = load_tx_model(gff)
    ref_gff = os.path.join(workdir, "01_genome", _seg(sp), "genomic.gff")
    gene_src = ref_gff if os.path.isfile(ref_gff) else gff
    _gtx, genes = load_tx_model(gene_src)
    log("  基因几何：%s（%d 基因）" % (os.path.basename(gene_src), len(genes)))
    genes_by_contig = collections.defaultdict(list)
    for gid, g in genes.items():
        g = dict(g)
        g["gene"] = gid
        genes_by_contig[g["seqid"]].append(g)

    hits = collections.defaultdict(list)
    for r in rows:
        pep = iln(r.get("peptide", ""))
        name = str(r.get("event_orf") or r.get("protein") or "").split("/")[0]
        pos0 = r.get("event_pos0")
        if not pep or pos0 is None or not (parse_euk_name(name) or _l1sf(name)):
            continue
        cls, gene = classify_peptide(name, pos0, len(pep), txs, genes_by_contig)
        gap = _min_gap(pep_blocks(name, pos0, len(pep), txs), genes_by_contig)
        hits[pep].append({"class": cls, "gene": gene, "tier": r.get("tier", ""),
                          "orf": name, "pos": pos0, "ref_gap": gap})
    main = normalize_peptide_level(hits)
    counts = summarize(main)
    log("  真核事件（肽级归一）：%s" % counts)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
        path = os.path.join(out_dir, "euk_events.tsv")
        with _open_w(workdir, path, "w", encoding="utf-8", newline="") as fh:
            fh.write("peptide\tclass\tgene\ttier\tn_events\tall_classes\torf\tpos\tref_gap\n")
            for pep in sorted(main):
                v = main[pep]
                fh.write("\t".join([
                    pep, v["class"], v.get("gene", ""), str(v.get("tier", "")),
                    str(v.get("n_events", 1)), "|".join(v.get("all_classes", [])),
                    v.get("orf", ""), str(v.get("pos", "")), str(v.get("ref_gap", "")),
                ]) + "\n")
        log("  事件表 → %s" % path)
    return main


# ---------------------------------------------------------------- P8 入口
def _spectra_dirs(root, sp):
    seg = _seg(sp)
    if not seg:
        return
    for sub in ("04_search_disc", "05b_research"):
        base = _within(root, sub, seg)
        if os.path.isdir(base):
            yield base


def attach_positions(root, rows, log=print):
    """给 KEEP 行补 ``event_orf`` / ``event_pos0``（只认 L0| / LJR| / L2| 三种层）。"""
    targets = {(r.get("species", ""), iln(r.get("peptide", ""))) for r in rows}
    found = {}
    for sp in sorted({s for s, _p in targets if s}):
        for base in _spectra_dirs(root, sp):
            for task in sorted(os.listdir(base)):
                seg = _seg(task)
                if not seg:
                    continue
                spec = _within(base, seg, "pFind-Filtered.spectra")
                if not os.path.isfile(spec):
                    continue
                for r in load_spectra(spec):
                    pep = iln(r.get("Sequence") or "")
                    if not pep or (sp, pep) not in targets:
                        continue
                    prots = [x.strip() for x in (r.get("Proteins") or "")
                             .replace(";", "/").split("/") if x.strip()]
                    poss = [x.strip() for x in (r.get("Positions") or "")
                            .replace(";", "/").split("/") if x.strip()]
                    for prot, pos in zip(prots, poss):
                        if not parse_euk_name(prot):
                            continue
                        try:
                            p0 = int(pos.split(",")[0])
                        except (ValueError, IndexError):
                            continue
                        found.setdefault((sp, pep), (prot, p0))
                        break
    n = 0
    for r in rows:
        key = (r.get("species", ""), iln(r.get("peptide", "")))
        hit = found.get(key)
        if hit:
            r["event_orf"], r["event_pos0"] = hit
            n += 1
        else:
            r.setdefault("event_orf", str(r.get("protein", "")).split("/")[0])
            r.setdefault("event_pos0", None)
    log("  真核坐标：%d/%d 行拿到层蛋白位点" % (n, len(rows)))
    return rows


def is_euk_project(root):
    """项目是否已用 eukbuild 建过真核层（02_db/<sp>/euk_manifest.json 存在）。"""
    db = Path(root, "02_db")
    if not db.is_dir():
        return False
    for p in db.iterdir():
        if p.is_dir() and (p / "euk_manifest.json").is_file():
            return True
    return False


def build_novgene_euk(root, out_dir, log=print):
    """P8 真核事件：KEEP 肽 → 七类 + 肽级归一，写 euk_events.tsv / summary。"""
    base = Path(out_dir).resolve()
    base.mkdir(parents=True, exist_ok=True)
    vdir = os.environ.get("ORFSCAPE_VALIDATE_DIR") or "06_validate"
    if "/" in vdir or "\\" in vdir or vdir in (".", ".."):
        vdir = "06_validate"
    src = Path(root, vdir, "novel_validated.tsv")
    rows_all = load_tsv(str(src)) if src.is_file() else []
    rows = [r for r in rows_all
            if str(r.get("final_keep", "")).startswith("KEEP")]
    log("  可信新肽 %d 行（过滤前 %d，不纳入统计）" % (len(rows), len(rows_all)))
    if not rows:
        log("  没有 KEEP 行，真核事件分类跳过")
        return {}

    attach_positions(root, rows, log=log)
    sps = sorted({r.get("species", "") for r in rows if r.get("species")})
    main = {}
    for sp in sps:
        sub = [r for r in rows if r.get("species") == sp]
        part = run(root, sp, sub, out_dir=None, log=log)
        for pep, v in part.items():
            v["species"] = sp
            main[(sp, pep)] = v

    counts = collections.Counter(v["class"] for v in main.values())
    genes = collections.defaultdict(set)
    for v in main.values():
        if v.get("gene"):
            genes[v["class"]].add(v["gene"])

    lines = ["species\tpeptide\tclass\tgene\ttier\tn_events\tall_classes\torf\tpos"
             "\tref_gap"]
    for (sp, pep) in sorted(main):
        v = main[(sp, pep)]
        lines.append("\t".join([
            sp, pep, v["class"], v.get("gene", ""), str(v.get("tier", "")),
            str(v.get("n_events", 1)), "|".join(v.get("all_classes", [])),
            v.get("orf", ""), str(v.get("pos", "")), str(v.get("ref_gap", "")),
        ]))
    (base / "euk_events.tsv").write_text("\n".join(lines) + "\n", encoding="utf-8")

    # 基因间/内含子没有可命名的基因符号，只统计基因数会出现 "0" 的误导；
    # 因此同时给 ORF 数（每个事件归属的 ORF 去重）
    orfs = collections.defaultdict(set)
    for v in main.values():
        if v.get("orf"):
            orfs[v["class"]].add(v["orf"])
    summ = ["class\tlabel_en\tpeptides\torfs\tgenes"]
    for cls in CLASSES:
        if not counts.get(cls):
            continue
        summ.append("%s\t%s\t%d\t%d\t%d" % (cls, CLASS_LABEL[cls][1],
                                            counts[cls], len(orfs.get(cls, ())),
                                            len(genes.get(cls, ()))))
    (base / "euk_events_summary.tsv").write_text("\n".join(summ) + "\n", encoding="utf-8")
    log("  真核事件（肽级归一）：%s" % dict(counts))
    log("  事件表 → %s" % (base / "euk_events.tsv"))
    return main
