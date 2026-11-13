# -*- coding: utf-8 -*-
"""新基因事件：按基因组坐标 + 同链判定分类（仅统计可信新肽）。

分类规则（用户口径 2026-11-06）：
  1. 新肽段先按 ORF 的坐标与链映射回基因组区间；
  2. 只与**同链**注释基因比较——同链重叠 → 已注释基因校正
     （肽段越过注释起点 = N 端延伸；否则 = 同链基因内读框/内部 ORF）；
  3. 无同链重叠 → 基因间区新编码基因；反义链重叠不改变分类
     （反义 ORF 是否新基因只看它自身所在链的注释状态）；
  4. 统计对象仅限**通过过滤的可信新肽**（final_keep=KEEP）。

坐标约定（基因组重构验证，含 I/L 归一）：
  L2_<contig>_<±k>_<start>：正链 residue1 = start-1+k；负链 residue1 = start-1-k
"""
import os
import re

import matplotlib.pyplot as plt
import numpy as np

from .novel import iln, l1_loci, load_spectra, read_faa
from .report import COLORS, COL1, COL2, _apply_zh_font, _save, load_tsv

TYPE_INTERGENIC = "intergenic_6frame"
TYPE_CORRECTED = "annotated_correction"
SUB_NTERM = "nterm_extension"
SUB_SAME_STRAND_FRAME = "same_strand_other_frame"
SUB_ANTISENSE_OVERLAP = "antisense_overlap"
SUB_INTERGENIC_ONLY = "intergenic_only"
SUB_INTERGENIC_OTHER = "unmapped"

TYPE_LABEL = {
    TYPE_INTERGENIC: ("基因间区新编码基因", "Intergenic new coding gene"),
    TYPE_CORRECTED: ("已注释基因校正", "Annotated gene correction"),
}
SUB_LABEL = {
    SUB_NTERM: ("N 端延伸（与注释 CDS 同读框）",
                "N-term extension (same frame as annotated CDS)"),
    SUB_SAME_STRAND_FRAME: ("新编码基因·注释基因区内不同读框",
                            "New gene · different frame inside a gene"),
    SUB_ANTISENSE_OVERLAP: ("新编码基因·反义链重叠",
                            "New gene · antisense overlap"),
    SUB_INTERGENIC_ONLY: ("新编码基因·纯基因间区",
                          "New gene · purely intergenic"),
    SUB_INTERGENIC_OTHER: ("新编码基因·未定位",
                           "New gene · unmapped"),
}


def T(zh, en):
    from .i18n import get_lang
    return zh if get_lang() == "zh" else en


def parse_orf_name(name):
    """L2_<contig>_<sign><off>_<start> → (contig, sign, off, start) 或 None。"""
    p = str(name or "").split("_")
    if len(p) < 4 or p[0] != "L2":
        return None
    frame = p[-2]
    if frame[0] not in "+-" or not frame[1:].isdigit():
        return None
    try:
        return "_".join(p[1:-2]), frame[0], int(frame[1:]), int(p[-1])
    except ValueError:
        return None


def parse_l6_name(name):
    """`L6|<contig>|<strand>|<start>-<end>|<n>` → (contig, strand, start1, end1)。

    转录证据层的坐标**已经是基因组坐标**（由 `dbbuild.l6_import` 回贴过），
    所以这里不需要再换算坐标系，直接用起点/终点算肽段跨度。
    """
    f = str(name or "").split("|")
    if len(f) < 4 or f[0] != "L6":
        return None
    if f[2] not in ("+", "-"):
        return None
    try:
        a, b = f[3].split("-", 1)
        return f[1], f[2], int(a), int(b)
    except ValueError:
        return None


def _species_list(rows):
    """从数据行推断物种（不写死 Mab/Mch）。"""
    sps = sorted({str(r.get("species", "") or "") for r in rows} - {""})
    return sps or ["Mab", "Mch"]


def parse_l3_name(name):
    """`L3|ext|<protein_id>|up<NN>aa` → (protein_id, upstream_aa) 或 None。"""
    f = str(name or "").split("|")
    if len(f) < 4 or f[0] != "L3":
        return None
    digits = "".join(ch for ch in f[3] if ch.isdigit())
    try:
        return f[2], int(digits)
    except ValueError:
        return None


def _peptide_span(name, pos1, pep_len, clen, l1_map=None):
    """肽段 → (contig, strand, lo, hi)（正链 0-based 闭区间）。L1/L2/L6 三种命名都支持。"""
    if str(name).startswith("L1_"):
        hit = (l1_map or {}).get(str(name))
        if not hit:
            return None
        contig, a, b, strand = hit
        off = int(pos1) * 3                # pFind Positions 是 0-based 氨基酸起点
        n3 = pep_len * 3
        if strand == "+":
            lo = a - 1 + off
            return contig, "+", lo, lo + n3 - 1
        hi = b - off
        return contig, "-", hi - n3 + 1, hi
    if str(name).startswith("L6|"):
        p = parse_l6_name(name)
        if not p:
            return None
        contig, strand, a, b = p
        off = int(pos1) * 3                # pFind Positions 是 0-based 氨基酸起点
        n3 = pep_len * 3
        if strand == "+":
            lo = a - 1 + off
            return contig, "+", lo, lo + n3 - 1
        hi = b - off
        return contig, "-", hi - n3 + 1, hi
    parsed = parse_orf_name(name)
    if not parsed or parsed[0] not in clen:
        return None
    lo, hi = peptide_plus_span(parsed, clen[parsed[0]], pos1, pep_len)
    return parsed[0], parsed[1], lo, hi


def _dirs(base):
    """base 下的物种子目录（不写死物种代码）。"""
    if not os.path.isdir(base):
        return []
    return sorted(x for x in os.listdir(base)
                  if os.path.isdir(os.path.join(base, x)))


def contig_lengths(root, sp):
    """contig 长度（genomic.fna）。"""
    fa = os.path.join(root, "01_genome", sp, "genomic.fna")
    out = {}
    if not os.path.exists(fa):
        return out
    name, n = None, 0
    with open(fa, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    out[name] = n
                name, n = line[1:].split()[0], 0
            else:
                n += len(line.strip())
    if name:
        out[name] = n
    return out


def orf_plus_span(parsed, clen, orf_aa_len):
    """ORF 在正链基因组上的 0-based 闭区间。

    正链：start 即首碱基；负链：start 是**反向互补坐标系**下的首碱基。
    """
    contig, sign, off, start = parsed
    n3 = orf_aa_len * 3
    if sign == "+":
        lo = start - 1
        return lo, lo + n3 - 1
    rc_lo = start - 1
    rc_hi = rc_lo + n3 - 1
    return clen - 1 - rc_hi, clen - 1 - rc_lo


def peptide_plus_span(parsed, clen, pos0, plen):
    """肽段（Positions 列 0-based 氨基酸起点）→ 正链 0-based 闭区间。"""
    contig, sign, off, start = parsed
    n3 = plen * 3
    if sign == "+":
        lo = start - 1 + pos0 * 3
        return lo, lo + n3 - 1
    rc_pos = start - 1 + pos0 * 3
    hi = clen - 1 - rc_pos
    return hi - n3 + 1, hi


def parse_gff(root, sp):
    """注释基因区间 → {contig: [(start1, end1, strand, locus_tag)]}。"""
    path = os.path.join(root, "01_genome", sp, "genomic.gff")
    out = {}
    if not os.path.exists(path):
        return out
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9 or c[2] != "gene":
                continue
            try:
                a, b = int(c[3]), int(c[4])
            except ValueError:
                continue
            m = re.search(r"locus_tag=([^;]+)", c[8])
            out.setdefault(c[0], []).append(
                (a, b, c[6], m.group(1) if m else ""))
    for k in out:
        out[k].sort()
    return out


def peptide_positions(root, targets):
    """{(species, peptide): [(L2 ORF, pos1)]} —— E2 + P6b 回搜谱表 Positions（1-based）。"""
    out = {}
    bases = [os.path.join(root, "04_search_disc")]
    # P6b 回搜目录（05b_research）必须**总是**扫：只在 P6b 里出现的肽
    # （例如只被小库回搜到的）否则拿不到位点，事件分类会误判"未定位"。
    # 以前靠 ORFSCAPE_VALIDATE_DIR 才扫，走其它入口时就会漏（TP 实测 14 行）。
    p6b = os.path.join(root, "05b_research")
    if os.path.isdir(p6b):
        bases.append(p6b)
    for base in bases:
        for sp in _dirs(base):
            disc = os.path.join(base, sp)
            if not os.path.isdir(disc):
                continue
            for task in sorted(os.listdir(disc)):
                spec = os.path.join(disc, task, "pFind-Filtered.spectra")
                if not os.path.exists(spec):
                    continue
                for r in load_spectra(spec):
                    pep = iln(r.get("Sequence") or "")
                    if not pep or (sp, pep) not in targets:
                        continue
                    prots = [x.strip() for x in (r.get("Proteins") or "")
                             .replace(";", "/").split("/") if x.strip()]
                    poss = [x.strip() for x in (r.get("Positions") or "")
                            .replace(";", "/").split("/") if x.strip()]
                    hit = None
                    for prot, pos in zip(prots, poss):
                        # L2 六框 / L6 转录 ORF 的位点可直接换算基因组坐标；
                        # L3 的位点是相对"延伸后蛋白"的，单独走 parse_l3_name 分支。
                        if not (prot.startswith("L2_") or prot.startswith("L6|")
                                or prot.startswith("L3")):
                            continue
                        try:
                            p1 = int(pos.split(",")[0])
                        except (ValueError, IndexError):
                            continue
                        hit = (prot, p1)
                        break            # Proteins 按打分排序，取第一个（最优）
                    if hit:
                        lst = out.setdefault((sp, pep), [])
                        if hit not in lst:
                            lst.append(hit)
    return out


def classify_events(root, rows, log=print):
    """可信新肽 → 事件类型。

    规则（用户口径）：肽段落在**与注释蛋白同读框**的六框 ORF 内 → N 端延伸
    （已注释基因校正）；其余（不同读框 / 反义链 / 纯基因间区）→ 新编码基因。
    """
    targets = {(r.get("species", ""), r.get("peptide", "")) for r in rows}
    sps = _species_list(rows)
    pos_map = peptide_positions(root, targets)
    cds = {sp: cds_features(root, sp) for sp in sps}
    clen = {sp: contig_lengths(root, sp) for sp in sps}
    l1_map = {}
    for sp in sps:
        l1_map.update(l1_loci(os.path.join(root, "02_db", sp)))
    out = []
    for r in rows:
        r = dict(r)
        sp, pep = r.get("species", ""), r.get("peptide", "")
        cand = pos_map.get((sp, pep), [])
        name, pos0 = (cand[0] if cand
                      else (str(r.get("protein", "")).split("/")[0], None))
        r["event_orf"] = name
        r["event_pos0"] = pos0
        etype, sub = TYPE_INTERGENIC, SUB_INTERGENIC_OTHER
        if str(name).startswith("L3"):
            # L3 是"注释 CDS + 上游延伸"：位点是相对延伸后蛋白的 1-based 氨基酸位。
            # 落在延伸段内 → 已注释基因校正（N 端延伸）；落在注释段内 → 同读框非同段。
            l3 = parse_l3_name(name)
            if l3 and pos0 is not None:
                _pid, up = l3
                etype = TYPE_CORRECTED
                # 延伸段占 0..up-1；pos0 是 0-based，起点 < up 才真正覆盖延伸段
                sub = SUB_NTERM if int(pos0) < up else SUB_SAME_STRAND_FRAME
        if pos0 is not None:
            span = _peptide_span(name, pos0, len(pep), clen.get(sp, {}), l1_map)
            if span:
                contig, strand, lo, hi = span
                same_cds = [c for c in cds.get(sp, {}).get(contig, [])
                            if c[2] == strand
                            and not (c[1] - 1 < lo or c[0] - 1 > hi)]
                anti = [c for c in cds.get(sp, {}).get(contig, [])
                        if c[2] != strand
                        and not (c[1] - 1 < lo or c[0] - 1 > hi)]
                if same_cds:
                    # 一个 ORF 常与多个注释 CDS 重叠（邻居、嵌套基因）。必须选**同框**的那个
                    # 来判定：旧代码取 same_cds[0]（按坐标第一个），会把"与 TPANIC_0425 同框、
                    # 向上延伸 41 aa"的肽错判成"不同框"（E005 就是这么错的）。
                    cand = []
                    for c in same_cds:
                        a, b, cstrand, ph = c
                        rel = ((lo - (a - 1) - ph) if cstrand == "+"
                               else ((b - 1) - ph - hi))
                        cand.append((rel, c))
                    inframe = [(r, c) for r, c in cand if r % 3 == 0]
                    if inframe:
                        # 同框：起点在注释起始上游（N 端延伸），或正好落在起始密码子上
                        # （非 ATG 起始被六框读成其氨基酸，同样说明注释起点未被使用）
                        etype, sub = TYPE_CORRECTED, SUB_NTERM
                    else:
                        sub = SUB_SAME_STRAND_FRAME
                elif anti:
                    sub = SUB_ANTISENSE_OVERLAP
                else:
                    sub = SUB_INTERGENIC_ONLY
        r["event_type"], r["event_sub"] = etype, sub
        out.append(r)
    kinds = {}
    for r in out:
        kinds[r["event_type"]] = kinds.get(r["event_type"], 0) + 1
    log(f"  坐标判定（可信肽）：{kinds}")
    return out


def cds_features(root, sp):
    """CDS 特征 → {contig: [(start1, end1, strand, phase)]}。"""
    path = os.path.join(root, "01_genome", sp, "genomic.gff")
    out = {}
    if not os.path.exists(path):
        return out
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9 or c[2] != "CDS":
                continue
            try:
                a, b, ph = int(c[3]), int(c[4]), int(c[7])
            except ValueError:
                continue
            out.setdefault(c[0], []).append((a, b, c[6], ph))
    for k in out:
        out[k].sort()
    return out


def peptide_samples(root):
    """肽 → 鉴定来源集合 {(species, sample_type)}（E2 的 L2/L3 命中 + P6b 回搜命中）。"""
    out = {}
    roots = [os.path.join(root, "04_search_disc")]
    # 同上：P6b 回搜目录总是扫，否则只在回搜里出现的肽会被算成"无样本归属"
    _p6b = os.path.join(root, "05b_research")
    if os.path.isdir(_p6b):
        roots.append(_p6b)
    for base in roots:
        if not os.path.isdir(base):
            continue
        for sp in _dirs(base):
            disc = os.path.join(base, sp)
            if not os.path.isdir(disc):
                continue
            for task in sorted(os.listdir(disc)):
                spec = os.path.join(disc, task, "pFind-Filtered.spectra")
                if not os.path.exists(spec):
                    continue
                parts = task.split("__")[0].split("_")
                if len(parts) < 3:
                    continue
                st = parts[1]
                for r in load_spectra(spec):
                    prot = (r.get("Proteins") or "").split(";")[0].rstrip("/")
                    if not (prot.startswith("L2_") or prot.startswith("L3")
                            or prot.startswith("L6|")):
                        continue
                    pep = iln(r.get("Sequence") or "")
                    if pep:
                        out.setdefault((sp, pep), set()).add((sp, st))
    return out


def summarize_events(rows, samples):
    """按 (样本, 事件类型/细分) 汇总肽段与蛋白数（样本归属非互斥）。"""
    def agg(sub):
        return dict(pep=len(sub),
                    prot=len({r.get("event_orf") or r.get("protein", "")
                              for r in sub}))
    res, subs = {}, {}
    for sp, st in samples:
        keep = [r for r in rows if (sp, st) in r.get("_samples", set())]
        res[(sp, st)] = {t: agg([r for r in keep if r["event_type"] == t])
                         for t in (TYPE_INTERGENIC, TYPE_CORRECTED)}
        subs[(sp, st)] = {k: agg([r for r in keep if r["event_sub"] == k])
                          for k in SUB_LABEL}
    res["ALL"] = {t: agg([r for r in rows if r["event_type"] == t])
                  for t in (TYPE_INTERGENIC, TYPE_CORRECTED)}
    subs["ALL"] = {k: agg([r for r in rows if r["event_sub"] == k])
                   for k in SUB_LABEL}
    return res, subs


def write_table(root, out_dir, counts, subs, samples):
    from .i18n import get_lang
    zh = get_lang() == "zh"
    cols = list(samples) + ["ALL"]
    head_cells = [T("事件种类", "Event type"), T("细分", "Subtype")]
    for c in cols:
        lab = T("合计", "Total") if c == "ALL" else c[0] + "·" + c[1]
        head_cells += [f"{lab} {T('肽段', 'pep')}", f"{lab} {T('蛋白', 'prot')}"]
    rows_out = [head_cells]
    subs_of = {TYPE_CORRECTED: [SUB_NTERM],
               TYPE_INTERGENIC: [SUB_SAME_STRAND_FRAME, SUB_ANTISENSE_OVERLAP,
                                 SUB_INTERGENIC_ONLY, SUB_INTERGENIC_OTHER]}
    for t, lab in TYPE_LABEL.items():
        row = [lab[0 if zh else 1], ""]
        for c in cols:
            d = counts[c][t]
            row += [f"{d['pep']:,}", f"{d['prot']:,}"]
        rows_out.append(row)
        for k in subs_of[t]:
            if not any(subs[c][k]["pep"] for c in cols):
                continue
            srow = ["  " + SUB_LABEL[k][0 if zh else 1], ""]
            for c in cols:
                d = subs[c][k]
                srow += [f"{d['pep']:,}", f"{d['prot']:,}"]
            rows_out.append(srow)
    tbl = os.path.join(out_dir, "NG_table.tsv")
    with open(tbl, "w", newline="", encoding="utf-8-sig") as fh:
        for r in rows_out:
            fh.write(chr(9).join(r) + chr(10))
    try:
        from openpyxl import Workbook
        from openpyxl.styles import Alignment, Font
        wb = Workbook()
        ws = wb.active
        ws.title = "events by sample"
        for r in rows_out:
            ws.append(r)
        for c in ws[1]:
            c.font = Font(name="Arial", size=9, bold=True)
        for row in ws.iter_rows(min_row=2):
            for c in row:
                c.font = Font(name="Arial", size=9)
            row[0].alignment = Alignment(horizontal="left")
            row[1].alignment = Alignment(horizontal="left")
        ws.column_dimensions["A"].width = 26
        ws.column_dimensions["B"].width = 32
        for i2 in range(3, len(rows_out[0]) + 1):
            ws.column_dimensions[chr(64 + i2)].width = 13
        book = os.path.join(out_dir, "NG_table.xlsx")
        wb.save(book)
    except ImportError:
        pass


def ng1_events_bar(root, out_dir, counts, samples):
    from .i18n import get_lang
    zh = get_lang() == "zh"
    cols = list(samples)
    x = np.arange(len(cols))
    fig, axes = plt.subplots(1, 2, figsize=(COL2, 2.5),
                             gridspec_kw={"wspace": 0.34})
    for ax, key, ttl in ((axes[0], "pep", T("肽段数", "Peptides")),
                         (axes[1], "prot", T("蛋白数", "Proteins"))):
        inter = [counts[c][TYPE_INTERGENIC][key] for c in cols]
        corr = [counts[c][TYPE_CORRECTED][key] for c in cols]
        ax.bar(x, inter, 0.55, color=COLORS[0],
               label=TYPE_LABEL[TYPE_INTERGENIC][0 if zh else 1])
        ax.bar(x, corr, 0.55, bottom=inter, color=COLORS[1],
               label=TYPE_LABEL[TYPE_CORRECTED][0 if zh else 1])
        top = max([a + b for a, b in zip(inter, corr)]) or 1
        for xi, a, b in zip(x, inter, corr):
            ax.text(xi, a + b + top * 0.02, f"{a + b:,}", ha="center",
                    fontsize=5)
            if b:
                ax.text(xi, a + b / 2, f"{b}", ha="center", va="center",
                        fontsize=4.6, color="white")
        ax.set_xticks(x)
        ax.set_xticklabels([c[0] + chr(10) + c[1] for c in cols], fontsize=5.4)
        ax.set_ylabel(ttl, fontsize=5.6)
        ax.set_title(ttl, fontsize=7, fontweight="bold", loc="left")
        ax.set_ylim(0, top * 1.32)
    axes[0].legend(frameon=False, fontsize=4.8, loc="upper right")
    fig.suptitle(T("新基因事件（按物种 × 鉴定来源）",
                   "Novel gene events by species x sample type"),
                 fontsize=7.4, fontweight="bold", x=0.01, ha="left", y=1.05)
    _save(fig, out_dir, "NG1_events_bar")


def peptide_nterm_mods(root, targets):
    """带肽段 N 端修饰（位置 0/含 N-term 标记）的肽集合 {(sp, pep)}。"""
    out = set()
    base = os.path.join(root, "04_search_disc")
    for sp in _dirs(base):
        disc = os.path.join(base, sp)
        if not os.path.isdir(disc):
            continue
        for task in sorted(os.listdir(disc)):
            spec = os.path.join(disc, task, "pFind-Filtered.spectra")
            if not os.path.exists(spec):
                continue
            for r in load_spectra(spec):
                pep = iln(r.get("Sequence") or "")
                if not pep or (sp, pep) not in targets:
                    continue
                mod = r.get("Modification") or ""
                for part in filter(None, mod.split(";")):
                    f = part.split(",")
                    if len(f) < 2:
                        continue
                    try:
                        pos = int(f[0])
                    except ValueError:
                        continue
                    if pos == 0 or "N-term" in f[1]:
                        out.add((sp, pep))
                        break
    return out


def nterm_evidence(sub_rows, orf_first, nmods, sp_list):
    """按三判据统计 ORF 的 N 端证据。

    判据：①肽段起点 = ORF 第 1 位；②起点 = 第 2 位且第 1 位为 M（起始 Met）；
          ③肽段带 N 端修饰。任一满足即该 ORF 有 N 端证据。
    """
    ev, detail = set(), {"pos1": set(), "met2": set(), "mod": set()}
    for r in sub_rows:
        orf = r.get("event_orf", "")
        sp = r.get("species", "")
        if sp not in sp_list:
            continue
        p0 = r.get("event_pos0")
        if p0 is None:
            continue
        p0 = int(p0)
        if p0 == 0:
            ev.add(orf)
            detail["pos1"].add(orf)
        elif p0 == 1 and orf_first.get(orf, "") == "M":
            ev.add(orf)
            detail["met2"].add(orf)
        if (sp, r.get("peptide", "")) in nmods:
            ev.add(orf)
            detail["mod"].add(orf)
    return ev, detail


def ng2_orf_distributions(root, out_dir, rows, lens, samples,
                          orf_first=None, nmods=None):
    """新编码基因分布：按样本分组的 ORF 长度 / 每基因肽段数 / N 端支持。"""
    cols = list(samples)
    sp_list = _species_list(rows)
    pal = [COLORS[0], COLORS[1], "#9CA3AF"]
    fig, axes = plt.subplots(1, 3, figsize=(COL2, 2.6),
                             gridspec_kw={"wspace": 0.40})
    ax = axes[0]
    data, labels = [], []
    for sp, st in cols:
        sub = [r for r in rows if r["event_type"] == TYPE_INTERGENIC
               and (sp, st) in r.get("_samples", set())]
        vals = [lens[r["event_orf"]] for r in sub if r["event_orf"] in lens]
        data.append(vals if vals else [np.nan])
        labels.append(sp + chr(10) + st)
    bp = ax.boxplot(data, labels=labels, widths=0.5, patch_artist=True,
                    showfliers=False,
                    medianprops=dict(color="#111827", lw=0.9))
    for i, patch in enumerate(bp["boxes"]):
        patch.set_facecolor(pal[i % len(pal)])
        patch.set_alpha(0.5)
        patch.set_edgecolor(pal[i % len(pal)])
    ax.set_ylabel(T("新 ORF 长度 (aa)", "ORF length (aa)"), fontsize=5.6)
    ax.set_title(T("a  新 ORF 长度分布", "a  ORF length distribution"),
                 fontsize=7, fontweight="bold", loc="left")
    ax = axes[1]
    cats = [(1, 1, "1"), (2, 2, "2"), (3, 5, "3-5"), (6, 999, ">=6")]
    w = 0.26
    for i, (sp, st) in enumerate(cols):
        sub = [r for r in rows if r["event_type"] == TYPE_INTERGENIC
               and (sp, st) in r.get("_samples", set())]
        per = {}
        for r in sub:
            k = r.get("event_orf", "")
            per[k] = per.get(k, 0) + 1
        vals = [sum(1 for v in per.values() if lo <= v <= hi)
                for lo, hi, _ in cats]
        xs = np.arange(len(cats)) + (i - 1) * w
        bars = ax.bar(xs, vals, w, color=pal[i % len(pal)], alpha=0.9)
        for b, v in zip(bars, vals):
            if v:
                ax.text(b.get_x() + b.get_width() / 2, v + 0.4, f"{v}",
                        ha="center", fontsize=4.4)
    ax.set_xticks(np.arange(len(cats)))
    ax.set_xticklabels([c[2] for c in cats], fontsize=6)
    ax.set_xlabel(T("每基因肽段数", "Peptides per gene"), fontsize=5.6)
    ax.set_ylabel(T("新编码基因数", "ORFs"), fontsize=5.6)
    ax.set_title(T("b  新编码基因肽段数分布", "b  Peptides per novel ORF"),
                 fontsize=7, fontweight="bold", loc="left")
    ax = axes[2]
    w = 0.26
    for i, (sp, st) in enumerate(cols):
        sub = [r for r in rows if r["event_type"] == TYPE_INTERGENIC
               and (sp, st) in r.get("_samples", set())]
        orfs = {r["event_orf"] for r in sub}
        ev, _ = nterm_evidence(sub, orf_first or {}, nmods or set(),
                               sp_list)
        nt = len(ev)
        vals = [nt, len(orfs) - nt]
        xs = np.arange(2) + (i - 1) * w
        bars = ax.bar(xs, vals, w, color=pal[i % len(pal)], alpha=0.9)
        for b, v in zip(bars, vals):
            ax.text(b.get_x() + b.get_width() / 2, v + 0.3, f"{v}",
                    ha="center", fontsize=4.4)
    ax.set_xticks(np.arange(2))
    ax.set_xticklabels([T("含 N 端肽段", "N-term"), T("仅内部", "Internal")],
                       fontsize=5.4)
    ax.set_ylabel(T("新编码基因数", "ORFs"), fontsize=5.6)
    ax.set_title(T("c  新编码基因 N 端预测支持", "c  N-term start support"),
                 fontsize=7, fontweight="bold", loc="left")
    handles = [plt.Rectangle((0, 0), 1, 1, fc=pal[i % len(pal)], alpha=0.9)
               for i in range(len(cols))]
    fig.legend(handles, [sp + "·" + st for sp, st in cols], frameon=False,
               fontsize=5, ncol=3, loc="upper right",
               bbox_to_anchor=(1.0, 1.07))
    _save(fig, out_dir, "NG2_orf_distributions")
    out = {}
    for sp, st in cols:
        sub = [r for r in rows if r["event_type"] == TYPE_INTERGENIC
               and (sp, st) in r.get("_samples", set())]
        orfs = {r["event_orf"] for r in sub}
        vals = [lens[r["event_orf"]] for r in sub if r["event_orf"] in lens]
        ev, _d = nterm_evidence(sub, orf_first or {}, nmods or set(),
                                sp_list)
        out[(sp, st)] = dict(
            n_orf=len(orfs),
            med_len=float(np.median(vals)) if vals else 0.0,
            nterm=len(ev))
    return out


def write_conclusions(root, out_dir, counts, subs, st_map, samples):
    from .i18n import get_lang
    zh = get_lang() == "zh"
    cols = list(samples)
    parts_zh = "；".join(
        f"{sp}·{st}: 新编码基因 {counts[(sp, st)][TYPE_INTERGENIC]['pep']:,} 肽 / "
        f"{counts[(sp, st)][TYPE_INTERGENIC]['prot']:,} ORF、N 端延伸 "
        f"{counts[(sp, st)][TYPE_CORRECTED]['pep']:,} 肽 / "
        f"{counts[(sp, st)][TYPE_CORRECTED]['prot']:,} 基因"
        for sp, st in cols)
    parts_en = "; ".join(
        f"{sp}·{st}: new genes "
        f"{counts[(sp, st)][TYPE_INTERGENIC]['pep']:,} pep / "
        f"{counts[(sp, st)][TYPE_INTERGENIC]['prot']:,} ORFs, extension "
        f"{counts[(sp, st)][TYPE_CORRECTED]['pep']:,} pep / "
        f"{counts[(sp, st)][TYPE_CORRECTED]['prot']:,} genes"
        for sp, st in cols)
    pa = counts["ALL"]
    if zh:
        lines = [
            "# 新基因事件 —— 结论句（物种 × 鉴定来源；仅可信新肽）",
            "",
            f"**NG1 · 事件分类（分物种/来源）** —— {parts_zh}；合计（肽级去重）"
            f"新编码基因 {pa[TYPE_INTERGENIC]['pep']:,} 肽 / "
            f"{pa[TYPE_INTERGENIC]['prot']:,} ORF、N 端延伸 "
            f"{pa[TYPE_CORRECTED]['pep']:,} 肽 / "
            f"{pa[TYPE_CORRECTED]['prot']:,} 基因。样本归属非互斥。",
            "",
            f"**NG1b · 新编码基因细分（合计）** —— " + "；".join(
                f"{SUB_LABEL[k][0]} {subs['ALL'][k]['pep']:,} 肽 / "
                f"{subs['ALL'][k]['prot']:,} ORF"
                for k in (SUB_SAME_STRAND_FRAME, SUB_ANTISENSE_OVERLAP,
                          SUB_INTERGENIC_ONLY)) + "。",
            "",
            f"**NG2 · 新编码基因分布（分样本）** —— N 端证据按三判据合并："
            f"肽段起点=ORF 第 1 位、或起点=第 2 位且第 1 位为 M（起始 Met 切除）、"
            f"或肽段带 N 端修饰（如 Acetyl[AnyN-term]）。结果：" + "；".join(
                f"{sp}·{st} {st_map[(sp, st)]['n_orf']} 个 ORF"
                f"（长度中位 {st_map[(sp, st)]['med_len']:.0f} aa，"
                f"N 端支持 {st_map[(sp, st)]['nterm']}）"
                for sp, st in cols) + "。",
            "",
        ]
    else:
        lines = [
            "# Novel gene events — conclusions (species x sample type; "
            "credible peptides only)",
            "",
            f"**NG1 · Event classes (per species/source)** — {parts_en}; "
            f"totals (peptide-level dedup) new genes "
            f"{pa[TYPE_INTERGENIC]['pep']:,} pep / "
            f"{pa[TYPE_INTERGENIC]['prot']:,} ORFs, extension "
            f"{pa[TYPE_CORRECTED]['pep']:,} pep / "
            f"{pa[TYPE_CORRECTED]['prot']:,} genes. Attribution is "
            f"non-exclusive.",
            "",
            f"**NG1b · New-gene subtypes (total)** — " + "; ".join(
                f"{SUB_LABEL[k][1]} {subs['ALL'][k]['pep']:,} pep / "
                f"{subs['ALL'][k]['prot']:,} ORFs"
                for k in (SUB_SAME_STRAND_FRAME, SUB_ANTISENSE_OVERLAP,
                          SUB_INTERGENIC_ONLY)) + ".",
            "",
            f"**NG2 · New coding genes (per sample)** — N-terminal evidence "
            f"combines three criteria: peptide starts at ORF residue 1, or at "
            f"residue 2 with residue 1 = M (initiator Met excision), or the "
            f"peptide carries an N-terminal modification (e.g. "
            f"Acetyl[AnyN-term]). " + "; ".join(
                f"{sp}·{st} {st_map[(sp, st)]['n_orf']} ORFs "
                f"(median {st_map[(sp, st)]['med_len']:.0f} aa, "
                f"N-term support {st_map[(sp, st)]['nterm']})"
                for sp, st in cols) + ".",
            "",
        ]
    out = os.path.join(out_dir, "NG_figure_conclusions.md")
    with open(out, "w", encoding="utf-8") as fh:
        fh.write(chr(10).join(lines))
    return out


def build_novgene(root, out_dir, log=print):
    _apply_zh_font()
    os.makedirs(out_dir, exist_ok=True)
    vdir = os.environ.get("ORFSCAPE_VALIDATE_DIR", "06_validate")
    src = os.path.join(root, vdir, "novel_validated.tsv")
    rows_all = load_tsv(src) if os.path.exists(src) else []
    rows = [r for r in rows_all
            if str(r.get("final_keep", "")).startswith("KEEP")]
    log(f"  可信新肽 {len(rows):,} 行（过滤前 {len(rows_all):,}，不纳入统计）")
    rows = classify_events(root, rows, log=log)
    pmap = peptide_samples(root)
    samples = sorted({k for v in pmap.values() for k in v})
    for r in rows:
        r["_samples"] = pmap.get((r.get("species", ""),
                                  r.get("peptide", "")), set())
    counts, subs = summarize_events(rows, samples)
    lens, orf_first = {}, {}
    for sp in _species_list(rows):
        # L2 六框（长度 + 首残基用于起始 Met 判据）
        faa = os.path.join(root, "02_db", sp, "L2_6frame.faa")
        if os.path.exists(faa):
            for h, seq in read_faa(faa):
                if h.startswith("L2_"):
                    name = h.split()[0]
                    lens[name] = len(seq)
                    orf_first[name] = seq[:1]
        # L6 转录 ORF 同样要进长度表，否则 NG2 的长度分布会把 L6 事件漏掉
        faa6 = os.path.join(root, "02_db", sp, "L6_tx_orfs_filt.faa")
        if os.path.exists(faa6):
            for h, seq in read_faa(faa6):
                name = h.split()[0]
                if name.startswith("L6|"):
                    lens[name] = len(seq)
                    orf_first[name] = seq[:1]
    nmods = peptide_nterm_mods(
        root, {(r.get("species", ""), r.get("peptide", "")) for r in rows})
    write_table(root, out_dir, counts, subs, samples)
    ng1_events_bar(root, out_dir, counts, samples)
    st_map = ng2_orf_distributions(root, out_dir, rows, lens, samples,
                                   orf_first, nmods)
    write_conclusions(root, out_dir, counts, subs, st_map, samples)
    pa = counts["ALL"]
    log(f"  事件合计：新编码基因 {pa[TYPE_INTERGENIC]['pep']:,} 肽 / "
        f"{pa[TYPE_INTERGENIC]['prot']:,} ORF；N 端延伸 "
        f"{pa[TYPE_CORRECTED]['pep']:,} 肽 / "
        f"{pa[TYPE_CORRECTED]['prot']:,} 基因")
    return counts, subs, st_map
