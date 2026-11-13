# -*- coding: utf-8 -*-
"""C 建库模块（设计书 v1.1/v1.2）：原核多源层级库 + 两段漏斗拆分库。

层：
  L0 注释 CDS（RefSeq faa，I/L 归一去重）
  L1 ab initio（Prodigal/Augustus/GeneMarkS 三源；缺工具时整层跳过并声明——检测器降级链）
  L2 六框 ORF（fna 六框 stop-to-stop ≥10aa；仅保留对 L0 消化集有 ≥2 新肽的 ORF）
  L3 替代起始 N 端延伸（GFF CDS 上游同框延伸；锚序列设计：注释 CDS + 上游延伸段）
  L4 酶自切（Trypsin/LysC/LysN/LysargiNase 常见自切肽，CON_* 对照）
  L5 污染（内置精简 cRAP 集）
拆分库（两段漏斗发现轮）：
  db_annot = L0+L4+L5；db_aug = +L3；db_6f = +L2；db_tx = +L6（本轮空）
  各库蛋白级反转 decoy（I/L 归一口径一致）。
新肽判定：双酶规则（TL=KR_C、LL=KR_N）下，I/L 归一后不在 L0 消化集的肽。
"""
import os
import re

# ---------------------------------------------------------------- 基础
CODON = {}
_B = "TCAG"
_AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
_i = 0
for a in _B:
    for b in _B:
        for c in _B:
            CODON[a + b + c] = _AA[_i]
            _i += 1

IL = str.maketrans("IL", "LL")  # I→L 归一


def il_norm(seq):
    return seq.upper().translate(IL)


def translate(nt):
    aa = []
    for i in range(0, len(nt) - 2, 3):
        aa.append(CODON.get(nt[i:i + 3], "X"))
    return "".join(aa)


def revcomp(nt):
    return nt.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def six_frame_orfs(nt, min_aa=10):
    """六框 stop-to-stop ORF → [(frame, start0, aa)]，含跨行首尾方向。"""
    out = []
    for strand, seq in (("+", nt), ("-", revcomp(nt))):
        for f in range(3):
            prot = translate(seq[f:])
            for m in re.finditer(r"[^*]+", prot):
                aa = m.group(0)
                if len(aa) >= min_aa:
                    out.append((f"{strand}{f}", f + 3 * m.start(), aa))
    return out


def digest(seq, rule="KR_C", min_len=7, max_len=50):
    """in silico 酶切（无漏切口径，新肽判定够用）。KR_C=K/R 后切（P 位阻断）；KR_N=K/R 前切。"""
    peptides = []
    if rule == "KR_C":
        cuts = [0]
        for i, ch in enumerate(seq[:-1]):
            if ch in "KR" and seq[i + 1] != "P":
                cuts.append(i + 1)
        cuts.append(len(seq))
        for a, b in zip(cuts, cuts[1:]):
            if min_len <= b - a <= max_len:
                peptides.append(seq[a:b])
    else:  # KR_N：K/R 前切（其后为 P 时仍切——pFind KR P N 口径：P 在酶切位点下游时阻断）
        cuts = [0]
        for i, ch in enumerate(seq):
            if i > 0 and ch in "KR":
                cuts.append(i)
        cuts.append(len(seq))
        for a, b in zip(cuts, cuts[1:]):
            if min_len <= b - a <= max_len:
                peptides.append(seq[a:b])
    return peptides


def read_faa(path):
    """[(header, seq)]。"""
    items = []
    h, buf = "", []
    for line in open(path, encoding="utf-8", errors="replace"):
        if line.startswith(">"):
            if h:
                items.append((h, "".join(buf)))
            h, buf = line[1:].strip(), []
        else:
            buf.append(line.strip())
    if h:
        items.append((h, "".join(buf)))
    return items


def write_faa(path, records):
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        for h, s in records:
            fh.write(">" + h + "\n")
            for i in range(0, len(s), 60):
                fh.write(s[i:i + 60] + "\n")


def read_fna(path):
    seqs, name, buf = {}, "", []
    for line in open(path, encoding="utf-8", errors="replace"):
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(buf)
            name = line[1:].strip().split()[0]
            buf = []
        else:
            buf.append(line.strip().upper())
    if name:
        seqs[name] = "".join(buf)
    return seqs


def parse_gff_cds(gff_path):
    """GFF3 → [dict(contig, start1, end1, strand, protein_id, locus_tag, partial)]。"""
    cds = []
    for line in open(gff_path, encoding="utf-8", errors="replace"):
        if line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        if len(p) < 9 or p[2] != "CDS":
            continue
        attrs = dict(kv.split("=", 1) for kv in p[8].split(";") if "=" in kv)
        if not attrs:                      # GTF 回退：key "value"; 形式
            for kv in p[8].split(";"):
                kv = kv.strip()
                parts = kv.split(None, 1)
                if len(parts) == 2:
                    attrs[parts[0]] = parts[1].strip().strip('"')
        if "pseudo" in attrs or attrs.get("exception", "") == "ribosomal_slippage":
            continue
        cds.append(dict(contig=p[0], start1=int(p[3]), end1=int(p[4]), strand=p[6] + "+-"[0],
                        protein_id=attrs.get("protein_id", ""),
                        locus_tag=attrs.get("locus_tag", attrs.get("Parent", "gene")[:-4] if attrs.get("Parent") else ""),
                        partial=("partial=true" in p[8] or "start_range" in p[8] or "end_range" in p[8])))
    for c in cds:
        c["strand"] = "+" if c["strand"][0] == "+" else "-"
    return cds


# ---------------------------------------------------------------- 酶自切与污染（内置精简集）
ENZYME_AUTOLYSIS = {
    "CON_Trypsin": ["LSSPATLNSR", "VATEESCDVR", "LGADSGQYGPFQ", "SGVSSSTGR", "IQSVESSR"],
    "CON_LysC": ["ADEEEEDDVLGK", "GEYHGGETAAELSSR", "TNEIVEALYK"],
    "CON_LysN": ["KPDGTFEEK", "KELGSDANK", "KLSIDEAR"],
    "CON_LysargiNase": ["RITVNVAEVIHELR", "RVLYVGSEAAR", "RALNTLIAK"],
}
CONTAMINANTS = {
    "Crap_KeratinK1": "MSTVSSQSSSSSSSVRSPSGGGGGGRGSSSGSSRSGSSSSSSSSSGG",
    "Crap_KeratinK10": "MSTSSSSSSSSSFVSSSSSSSSSGGSSSSSSSSSSGRSGGSSSSSSSGG",
    "Crap_BSA": "MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPFDEHVKLVNELTEFAKTCVADESHAGCEKSLHTLFGDELCKVASLRETYGDMADCCEKQEPERNECFLSHKDDSPDLPKLKPDPNTLCDEFKADEKKFWGKYLYEIARRHPYFYAPELLYYANKYNGVFQECCQAEDKGACLLPKIETMREKVLASSARQRLRCASIQKFGERALKAWSVARLSQKFPKAEFVEVTKLVTDLTKVHKECCHGDLLECADDRADLAKYICDNQDTISSKLKECCDKPLLEKSHCIAEVEKDAIPENLPPLTADFAEDKDVCKNYQEAKDAFLGSFLYEYSRRHPEYAVSVLLRLAKEYEATLEECCAKDDPHACYSTVFDKLKHLVDEPQNLIKQNCDQFEKLGEYGFQNALIVRYTRKVPQVSTPTLVEVSRSLGKVGTRCCTKPESERMPCTEDYLSLILNRLCVLHEKTPVSDRVTKCCTESLVNRRPCFSALEVDETYVPKEFNAETFTFHADICTLSEKERQIKKQTALVELVKHKPKATKEQLKAVMDDFAAFVEKCCKADDKETCFAEEGKKLVAASQAALGL",
}


# ---------------------------------------------------------------- L1 / L6
def _find_prodigal():
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    cands = [os.path.join(here, "tools", "prodigal.exe"),
             os.path.join(here, "tools", "prodigal")]
    for c in cands:
        if os.path.exists(c):
            return c
    import shutil
    return shutil.which("prodigal")


def _overlap(a0, a1, b0, b1):
    lo, hi = max(a0, b0), min(a1, b1)
    return max(0, hi - lo + 1)


def build_l1(fna, out_dir, cds_list, log=print, min_aa=30, max_cds_overlap=0.8):
    """Prodigal ab initio 层。返回 (filt_records, note)。

    过滤：< min_aa 丢弃；与同链注释 CDS 重叠超过 max_cds_overlap × ORF 长度丢弃
    （这些是注释基因的等价物，留在 L1 只会重复计分）。
    """
    exe = _find_prodigal()
    if not exe:
        return [], "L1 跳过（未找到 Prodigal）"
    faa = os.path.join(out_dir, "L1_prodigal.faa")
    gff = os.path.join(out_dir, "L1_prodigal.gff")
    try:
        import subprocess
        cmd = [exe, "-i", fna, "-a", faa, "-f", "gff", "-o", gff,
               "-p", "single", "-g", "11"]
        subprocess.run(cmd, check=True, capture_output=True, timeout=1800)
    except Exception as e:
        return [], "L1 跳过（Prodigal 运行失败：%s）" % e

    # 坐标（prodigal gff：seqid / tool / CDS / start / end / . / strand / . / id=..;partial=..）
    coords = []
    try:
        for line in open(gff, encoding="utf-8", errors="replace"):
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[2] != "CDS":
                continue
            coords.append((p[0], int(p[3]), int(p[4]), p[6]))
    except OSError:
        coords = []

    ann = [(c["start1"], c["end1"], c["strand"]) for c in cds_list] if cds_list else []
    recs = read_faa(faa)
    kept = []
    keep_coords = []
    for i, (h, s) in enumerate(recs):
        if len(s) < min_aa:
            continue
        if ann and i < len(coords):
            ctg, s0, e0, st = coords[i]
            if any(st == ast and _overlap(s0, e0, as0, as1) >
                   max_cds_overlap * (e0 - s0 + 1) for as0, as1, ast in ann):
                continue
        name = f"L1_prodigal_{i + 1}"
        kept.append((name, il_norm(s)))
        # L1 的坐标必须落盘：P6 的 locus 归并与 P8 的事件判定都要它，
        # 只留名字的话这一层只能当"无位置的新肽"，事件分类会全判未定位。
        if i < len(coords):
            ctg, s0, e0, st = coords[i]
            keep_coords.append((name, ctg, s0, e0, st))
    write_faa(os.path.join(out_dir, "L1_prodigal_filt.faa"), kept)
    _write_l1_coords(out_dir, keep_coords, append=False)
    return kept, "L1 Prodigal：%d 条 → 过滤后 %d 条" % (len(recs), len(kept))


def _write_l1_coords(out_dir, rows, append=True):
    """L1 坐标表：`<orf>\\t<contig>\\t<start>\\t<end>\\t<strand>`（1-based 闭区间）。"""
    path = os.path.join(out_dir, "L1_coords.tsv")
    mode = "a" if append and os.path.exists(path) else "w"
    with open(path, mode, encoding="utf-8", newline="\n") as fh:
        if mode == "w":
            fh.write("orf\tcontig\tstart\tend\tstrand\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")


def split_by_length(records, k=3):
    """按序列长度把记录**近似等分**成 k 池（按条数等分，不按阈值切）。

    边界取长序的 i/k 分位；长度恰等于边界时归入下一池（`>` 判据），
    因此各池条数可能相差百分之几，但不会出现 92%/8% 那种失衡。
    """
    if not records:
        return [[] for _ in range(k)]
    lens = sorted(len(s) for _, s in records)
    n = len(lens)
    bounds = [lens[int(n * i / k)] for i in range(1, k)]
    pools = [[] for _ in range(k)]
    for h, s in records:
        idx = sum(1 for b in bounds if len(s) > b)
        pools[min(idx, k - 1)].append((h, s))
    return pools


L1_EXTRA = ("L1_gms2_filt.faa", "L1_augustus_filt.faa", "L1_gms_filt.faa")


def load_extra_l1(out_dir, log=print):
    """读外部 ab initio 源（GeneMarkS-2 / Augustus）的过滤结果，按序列去重。

    这些文件由外部工具产出后放到 `02_db/<sp>/`（命名见 L1_EXTRA），
    本函数只做读取与去重，不调用任何 Linux 工具。
    """
    seen = set()
    out, notes = [], []
    for name in L1_EXTRA:
        p = os.path.join(out_dir, name)
        if not os.path.exists(p):
            continue
        recs = read_faa(p)
        n0 = len(recs)
        added = 0
        coords = []
        for h, s in recs:
            key = il_norm(s)
            if key in seen:
                continue
            seen.add(key)
            out.append((h, key))
            added += 1
            # 外部源若在头里带坐标（`>L1_gms2_1 CP004010.2 8341 9282 - ...`），
            # 一并收进 L1_coords.tsv，P6/P8 才能给 L1 命中定位。
            f = h.split()
            if len(f) >= 5 and f[4] in ("+", "-"):
                try:
                    coords.append((f[0], f[1], int(f[2]), int(f[3]), f[4]))
                except ValueError:
                    pass
        if coords:
            _write_l1_coords(out_dir, coords, append=True)
        notes.append("%s %d 条（去重后新入 %d）" % (name, n0, added))
        log("  " + notes[-1])
    if not out:
        notes.append("无外部 ab initio 源")
    return out, "；".join(notes)


def build_l6(out_dir, log=print, min_aa=30):
    """转录证据层：读 L6_tx_orfs.faa（由 RNA-seq 流程产出），过滤后返回。

    蛋白名**必须原样保留** `L6|<contig>|<strand>|<start>-<end>|_<n>`：
    P6（novel.py）用层白名单 `startswith("L6|")` 过滤，并用这段解析 locus
    `tx:<contig>:<floor200>`。改名成 `L6_tx_N` 会让全部 L6 命中静默丢弃。
    """
    src = os.path.join(out_dir, "L6_tx_orfs.faa")
    if not os.path.exists(src):
        return [], "L6 跳过（无 L6_tx_orfs.faa）"
    recs = read_faa(src)
    kept, bad, seen = [], 0, {}
    for h, s in recs:
        if len(s) < min_aa:
            continue
        rid = h.split()[0]
        if not rid.startswith("L6|"):
            bad += 1
            rid = "L6|unknown|+|0-0|%d" % bad   # 兜底：保住前缀与 4 段结构
        f = rid.split("|")
        if len(f) < 4 or "-" not in f[3]:
            bad += 1
            rid = "L6|unknown|+|0-0|%d" % bad
        n = seen.get(rid, 0) + 1                # Prodigal 同区域多 ORF 可能重名
        seen[rid] = n
        if n > 1:
            rid = "%s.%d" % (rid, n)
        kept.append((rid, il_norm(s)))
    write_faa(os.path.join(out_dir, "L6_tx_orfs_filt.faa"), kept)
    note = "L6 转录 ORF：%d → %d" % (len(recs), len(kept))
    if bad:
        note += "（%d 条头不合规已兜底）" % bad
    return kept, note


# ---------------------------------------------------------------- L6 坐标回贴
# Prodigal 是在"未注释转录区段序列"上跑的，它给的是区段内坐标（1-based）。
# 区段名 `L6|<contig>|<strand>|<s>-<e>|()` 只记了转录区本身，而取序列时
# 两端各扩了 50 bp，且 bedtools getfasta -s 对负链会反向互补 → 必须按链向换算。
def l6_parse_region(name):
    f = name.split("|")
    if len(f) < 4 or "-" not in f[3]:
        return None
    a, b = f[3].split("-", 1)
    try:
        return f[1], f[2], int(a), int(b)
    except ValueError:
        return None


def l6_map_cds(gff):
    """注释 GFF → {contig: [(start, end, strand)]}（只取 CDS，1-based 闭区间）。"""
    out = {}
    for line in open(gff, encoding="utf-8", errors="replace"):
        if line.startswith("#"):
            continue
        c = line.rstrip("\n").split("\t")
        if len(c) < 8 or c[2] != "CDS":
            continue
        try:
            out.setdefault(c[0], []).append((int(c[3]), int(c[4]), c[6]))
        except ValueError:
            continue
    return out


def l6_import(gff, faa, cds_gff, out_dir, min_aa=30, max_cds_overlap=0.8,
              log=print):
    """Prodigal 区段内坐标 → 基因组坐标；剔 <min_aa 与同链重叠 >阈值 的 ORF。

    产物 `L6_tx_orfs.faa`，头为 `L6|<contig>|<strand>|<start>-<end>|<n>`，
    与 novel.py 的层白名单 / locus 解析一致。
    """
    recs = read_faa(faa)
    cds = l6_map_cds(cds_gff)
    rows, bad = [], 0
    for h, seq in recs:
        head, _, tail = h.partition(" # ")
        nums = [x.strip() for x in ("# " + tail).split("#")[1:-1]] if tail else []
        reg = l6_parse_region(head)
        if len(nums) < 3 or reg is None:
            bad += 1
            continue
        try:
            p1, p2, pstrand = int(nums[0]), int(nums[1]), int(nums[2])
        except ValueError:
            bad += 1
            continue
        contig, rstrand, rs, re_ = reg
        g0 = max(0, rs - 1 - 50)
        g1 = re_ + 50
        if not (1 <= p1 <= p2 <= g1 - g0):
            bad += 1
            continue
        # 区段内区间 [p1,p2] → 基因组 1-based 闭区间；链向由"区段链"与
        # Prodigal 报的基因链共同决定（meta 模式在正链区段里也会报负链基因）。
        if rstrand == "-":
            gs, ge = g1 - p2 + 1, g1 - p1 + 1
        else:
            gs, ge = g0 + p1, g0 + p2
        strand = rstrand if pstrand > 0 else ("+" if rstrand == "-" else "-")
        rows.append((contig, gs, ge, strand, seq.upper().rstrip("*")))

    kept, n_short, n_stop, n_ovl = [], 0, 0, 0
    for contig, gs, ge, strand, seq in rows:
        if len(seq) < min_aa:
            n_short += 1
            continue
        if "*" in seq.rstrip("*"):
            n_stop += 1
            continue
        ovl = 0.0
        for s, e, st in cds.get(contig, []):
            if st != strand:
                continue
            lo, hi = max(gs, s), min(ge, e)
            if hi >= lo:
                ovl = max(ovl, (hi - lo + 1) / float(ge - gs + 1))
        if ovl > max_cds_overlap:
            n_ovl += 1
            continue
        kept.append((contig, gs, ge, strand, seq))

    kept.sort(key=lambda x: (x[0], x[1]))
    named = [("L6|%s|%s|%d-%d|%d" % (c, st, gs, ge, i), s)
             for i, (c, gs, ge, st, s) in enumerate(kept, 1)]
    write_faa(os.path.join(out_dir, "L6_tx_orfs.faa"), named)
    log("  L6 区段内 ORF %d → 入层 %d（<%d aa %d / 内部终止 %d / 同链 CDS>%d%% %d）"
        % (len(rows), len(named), min_aa, n_short, n_stop,
           int(max_cds_overlap * 100), n_ovl))
    if bad:
        log("  ！%d 条 Prodigal 头无法解析，已跳过" % bad)
    return named


# ---------------------------------------------------------------- 主构建
def build_species(code, genome_dir, out_dir, min_orf_aa=10, log=print):
    """对一个物种构建 L0–L5 + 拆分库。返回统计 dict。"""
    os.makedirs(out_dir, exist_ok=True)
    faa = os.path.join(genome_dir, "protein.faa")
    fna = os.path.join(genome_dir, "genomic.fna")
    gff = os.path.join(genome_dir, "genomic.gff")

    # L0：I/L 归一去重
    recs = read_faa(faa)
    l0, seen = [], set()
    for h, s in recs:
        s2 = il_norm(s)
        if s2 in seen:
            continue
        seen.add(s2)
        acc = h.split()[0]
        l0.append((f"L0_{acc}", s2))
    # L0 双酶消化集（新肽排除器）
    l0_peps = set()
    for _, s in l0:
        for rule in ("KR_C", "KR_N"):
            for p in digest(s, rule):
                l0_peps.add(p)

    # L2：六框 ORF（保 ≥2 新肽）
    genome = read_fna(fna)
    l2 = []
    for contig, nt in genome.items():
        for frame, start, aa in six_frame_orfs(nt, min_orf_aa):
            aail = il_norm(aa)
            newp = 0
            for rule in ("KR_C", "KR_N"):
                newp += sum(1 for p in digest(aail, rule) if p not in l0_peps)
                if newp >= 2:
                    break
            if newp >= 2:
                l2.append((f"L2_{contig}_{frame}_{start + 1}", aail))

    # L3：替代起始 N 端延伸（注释 CDS 上游同框，≥5 aa 且带新肽）
    l3 = []
    try:
        cds_list = parse_gff_cds(gff)
    except Exception as e:
        cds_list = []
        log(f"  ⚠ GFF 解析失败（L3 跳过）: {e}")
    for c in cds_list[:]:
        if c["partial"]:
            continue
        contig = c["contig"]
        nt = genome.get(contig)
        if not nt:
            continue
        s0, e0, strand = c["start1"] - 1, c["end1"], c["strand"]
        # 上游序列（按链方向），同框步进 3
        if strand == "+":
            up = nt[:s0]
            frame_off = (s0 - len(up))  # 0
            up_stream = nt[max(0, s0 - 3000):s0]
        else:
            up_stream = revcomp(nt[e0:e0 + 3000])
        aa_up = translate(up_stream[-(len(up_stream) // 3 * 3):])
        # 找最近的上游 stop，取其后的延伸段
        seg = aa_up[::-1]
        m = re.search(r"\*", seg)
        ext = seg[:m.start()] if m else seg[:50]
        ext = ext[::-1]
        if len(ext) < 5:
            continue
        extil = il_norm(ext)
        # 锚序列 = 延伸段 + 注释蛋白前 30aa（近似同框拼接；注释蛋白本身在 L0）
        anchor = l0_prot_prefix(c["protein_id"], l0, 30)
        if not anchor:
            continue
        extfull = extil + il_norm(anchor)
        has_new = any(p not in l0_peps for p in digest(extil, "KR_N")) or \
                  any(p not in l0_peps for p in digest(extil, "KR_C"))
        if has_new:
            l3.append((f"L3|ext|{c['protein_id']}|up{len(ext)}aa", extfull))

    # L1 ab initio（Prodigal；与注释 CDS 高度重叠的丢弃）
    l1, note_l1 = build_l1(fna, out_dir, cds_list, log=log)
    log("  " + note_l1)
    # L1 外部源（GeneMarkS-2 / Augustus 的过滤结果，放 02_db/<sp>/ 下即自动并入）
    l1_extra, note_extra = load_extra_l1(out_dir, log=log)
    if l1_extra:
        have = {il_norm(s) for _, s in l1}
        add = [(h, s) for h, s in l1_extra if s not in have]
        l1 = l1 + add
        note_l1 = "%s；外部源并入 %d（%s）" % (note_l1, len(add), note_extra)

    # L6 转录证据（若已由 RNA-seq 流程产出 L6_tx_orfs.faa）
    l6, note_l6 = build_l6(out_dir, log=log)
    log("  " + note_l6)

    # L4/L5
    l4 = [(f"{k}_autolysis{i}", il_norm(p)) for k, ps in ENZYME_AUTOLYSIS.items()
          for i, p in enumerate(ps, 1)]
    l5 = [(k, il_norm(s)) for k, s in CONTAMINANTS.items()]

    # 落盘各层
    write_faa(os.path.join(out_dir, "L0_annot.faa"), l0)
    write_faa(os.path.join(out_dir, "L2_6frame.faa"), l2)
    write_faa(os.path.join(out_dir, "L3_altstart.faa"), l3)
    write_faa(os.path.join(out_dir, "L4_autolysis.faa"), l4)
    write_faa(os.path.join(out_dir, "L5_contaminant.faa"), l5)

    # 拆分库（发现轮）+ 蛋白级反转 decoy
    def merged(name, layers):
        recs2 = [r for layer in layers for r in layer]
        write_faa(os.path.join(out_dir, name), recs2)
        dec = [("REV_" + h.replace(" ", "_")[:80], s[::-1]) for h, s in recs2]
        write_faa(os.path.join(out_dir, name.replace(".faa", "_decoy.faa")), dec)
        return len(recs2)
    base = [l0, l4, l5, l3, l1, l6]
    n_annot = merged("db_annot.faa", [l0, l4, l5])
    n_aug = merged("db_aug.faa", base)
    n_6f = merged("db_6f.faa", base + [l2])
    n_tx = merged("db_tx.faa", base)
    # 六框按长度**近似等分**成 k 个池（避免长 ORF 在竞争里压死短 ORF）
    pools = split_by_length(l2, k=3)
    n_pools = [merged("db_6f%d.faa" % (i + 1), base + [p])
               for i, p in enumerate(pools)]
    stats = dict(code=code, L0=len(l0), L1=len(l1), L2=len(l2), L3=len(l3),
                 L4=len(l4), L5=len(l5), L6=len(l6),
                 L2_pool1=len(pools[0]), L2_pool2=len(pools[1]),
                 L2_pool3=len(pools[2]),
                 db_annot=n_annot, db_aug=n_aug, db_6f=n_6f, db_tx=n_tx,
                 db_6f1=n_pools[0], db_6f2=n_pools[1], db_6f3=n_pools[2],
                 l0_peptides=len(l0_peps),
                 notes="%s；%s；六框按长度等分 3 池" % (note_l1, note_l6))
    return stats


def l0_prot_prefix(protein_id, l0, n):
    for h, s in l0:
        if protein_id and protein_id in h:
            return s[:n]
    return ""
