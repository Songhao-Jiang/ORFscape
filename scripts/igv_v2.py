# -*- coding: utf-8 -*-
"""IGV 风格事件页（v2，按 pdf_toolkit/README_case_pdf_mode.md 的版式策略）。

版式（自上而下，数据坐标 y 越大越靠上）：
    CDS 行（可多行）      y = base_y + i
    novel ORF 条          y = base_y - 1.0
    镜像核心浅阴影        跨 ORF 下方两条肽段行
    Trypsin 肽段行        y = base_y - 1.9
    LysargiNase 肽段行    y = base_y - 2.5
两面板：上 overview（ORF ± 25%，≥1.5 kb）/ 下 zoom（ORF + 6%，≤6 kb）。
"""
import math
import re

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.markers import CARETLEFT, CARETRIGHT


def _pack(items, gap):
    rows = []
    for it in sorted(items, key=lambda c: c["s"]):
        for row in rows:
            if it["s"] > row[-1]["e"] + gap:
                row.append(it)
                break
        else:
            rows.append([it])
    return rows


def _marker(strand):
    return CARETRIGHT if strand == "+" else CARETLEFT


def _chevrons(ax, a, b, y, strand, ax_w_pt, step_pt=15.0, size=2.6,
              color="black", alpha=1.0, lw=0.6):
    if b <= a:
        return
    x0, x1 = ax.get_xlim()
    step = max(step_pt / ax_w_pt * (x1 - x0), (x1 - x0) * 1e-6)
    # 箭头本身有宽度：横向内缩半个箭头，保证不越出框的边界
    inset = (size * 0.95) / max(ax_w_pt, 1.0) * (x1 - x0)
    lo, hi = a + inset, b - inset
    if hi <= lo:
        lo, hi = a, b
    n = max(1, int((hi - lo) / step))
    xs = [lo + (hi - lo) * (k + 0.5) / n for k in range(n)]
    ax.plot(xs, [y] * n, linestyle="none", marker=_marker(strand),
            markersize=size, markeredgewidth=lw, color=color, alpha=alpha,
            clip_on=True, zorder=4)


def _int_ticks(ax, x0, x1):
    span = max(x1 - x0, 1)
    step = 10 ** int(math.floor(math.log10(max(span / 6.0, 1))))
    for mult in (1, 2, 5, 10):
        if span / (step * mult) <= 7:
            step *= mult
            break
    ticks = [t for t in range(int(x0 // step) * step, int(x1) + step, step)
             if x0 <= t <= x1]
    ax.set_xticks(ticks)
    ax.set_xticklabels(["{:,}".format(t) for t in ticks])


def _frame_on_contig(start, end, strand, clen):
    """读框号（与 L2 命名同基准：正链从 contig 首碱基数，负链从反向互补首碱基数）。"""
    if strand == "+":
        return 0 if start is None else (int(start) - 1) % 3
    if not clen:
        return None
    return (int(clen) - int(end)) % 3


def _frame_txt(start, end, strand, clen, fallback=None):
    f = _frame_on_contig(start, end, strand, clen)
    if f is None:
        f = fallback
    return "" if f is None else "frame %s%d" % (strand, f)


def _draw_cds(ax, cds, x0, x1, base_y, ax_w_pt, orf_s, orf_strand,
              show_frame=True, label_fs=5.6, clen=None):
    sel = [c for c in cds if c["e"] >= x0 and c["s"] <= x1]
    if not sel:
        return 0, []
    rows = _pack(sel, gap=max(30, int((x1 - x0) * 0.004)))
    shapes = []
    for i, row in enumerate(rows):
        y = base_y + i
        for c in row:
            w = c["e"] - c["s"] + 1
            ax.add_patch(Rectangle((c["s"], y - 0.18), w, 0.36,
                                   facecolor="#E4E4E4", edgecolor="black",
                                   lw=0.6, clip_on=True, zorder=2))
            shapes.append((c["s"], c["e"], y, c["strand"]))
            lab = (c["locus"] + ("  " + c["gene"] if c["gene"] else "")).strip()
            frame_s = ""
            if show_frame and c["strand"] == orf_strand:
                frame_s = _frame_txt(c["s"], c["e"], c["strand"], clen,
                                     fallback=(c["s"] - orf_s) % 3)
            cx = min(max(c["s"] + w / 2.0, x0 + (x1 - x0) * 0.03),
                     x1 - (x1 - x0) * 0.03)
            # 用"点偏移"定位标签（不受 ylim/裁剪影响）：基因名离开条 4 pt，功能再上 5 pt
            ax.annotate("%s (%s)%s" % (lab, c["strand"],
                                       ("  " + frame_s) if frame_s else ""),
                        xy=(cx, y + 0.18), xycoords="data",
                        xytext=(0, 8), textcoords="offset points",
                        ha="center", va="bottom", fontsize=label_fs,
                        color="#1F2A38", annotation_clip=False, zorder=6)
            if c["product"]:
                ax.annotate(c["product"][:30],
                            xy=(cx, y + 0.18), xycoords="data",
                            xytext=(0, 3), textcoords="offset points",
                            ha="center", va="bottom", fontsize=4.5,
                            style="italic", color="#5A6B85",
                            annotation_clip=False, zorder=6)
    return len(rows), shapes


def _orf_label(frame_txt):
    """左栏标签：两行，避免超出左边界。"""
    return "novel ORF\n" + frame_txt if frame_txt else "novel ORF"


def _arrows(ax, shapes, ax_w_pt, ylo, yhi, ax_h_pt, step_pt=18.0, want=5.6):
    """按"框高（pt）× 0.8"定箭头大小：比原来 2.3/2.6 大一倍多，又不越出框。"""
    band = ax_h_pt * 0.36 / max(1e-9, (yhi - ylo))
    size = max(2.2, min(want, band * 0.72))
    for a, b, y, st in shapes:
        _chevrons(ax, a, b, y, st, ax_w_pt, step_pt=step_pt, size=size)
    return size


def _left_labels(ax, items, ylim):
    """左栏文字：把数据 y 映射成 axes 分数。"""
    lo, hi = ylim
    for txt, y, color in items:
        frac = (y - lo) / (hi - lo)
        ax.text(-0.012, frac, txt, transform=ax.transAxes, ha="right",
                va="center", fontsize=6.2, color=color)


def _mean_depth(segs, s, e):
    """区间平均深度（由 bedgraph 段加权）。"""
    if not segs or e < s:
        return None
    tot = 0.0
    for a, b, d in segs:
        lo, hi = max(a, s), min(b, e)
        if hi >= lo:
            tot += d * (hi - lo + 1)
    return tot / max(1, e - s + 1)


def _draw_depth(ax, segs, x0, x1, y_base, h=0.50, color="#5B8DB8"):
    """真·逐碱基深度（bedtools genomecov -bga 产出的 bedgraph）。"""
    win = [s for s in (segs or []) if s[1] >= x0 and s[0] <= x1]
    if not win or max(s[2] for s in win) <= 0:
        return 0.0
    top = max(s[2] for s in win)
    xs, ys = [], []
    for s0, e0, d in sorted(win):
        a, b = max(s0, x0), min(e0, x1)
        if b < a:
            continue
        xs += [a, b]
        ys += [y_base + h * d / top] * 2
    ax.fill_between(xs, y_base, ys, step="post", color=color, alpha=0.45,
                    linewidth=0, zorder=1)
    ax.plot(xs, ys, drawstyle="steps-post", color=color, lw=0.6, zorder=2)
    ax.annotate("max %.0f×" % top, xy=(x1, y_base + h * 0.98), xycoords="data",
                xytext=(-2, 1), textcoords="offset points", ha="right",
                va="bottom", fontsize=4.3, color="#2C5E86",
                annotation_clip=False, zorder=6)
    return top


def _draw_coverage(ax, models, contig, x0, x1, y_base, h=0.55, color="#5B8DB8"):
    """IGV 风格覆盖度：窗口内逐碱基深度（StringTie 模型深度）填充成阶梯面积。"""
    from orfscape.core.rnaabund import coverage_profile
    xs, ys = coverage_profile(models or [], contig, x0, x1)
    if not ys or max(ys) <= 0:
        return 0.0
    top = max(ys)
    yy = [y_base + h * v / top for v in ys]
    ax.fill_between(xs, y_base, yy, step="post", color=color, alpha=0.45,
                    linewidth=0, zorder=1)
    ax.plot(xs, yy, drawstyle="steps-post", color=color, lw=0.6, zorder=2)
    ax.annotate("max %.1f× (model)" % top, xy=(x1, y_base + h * 0.99), xycoords="data",
                xytext=(-2, 1), textcoords="offset points", ha="right",
                va="bottom", fontsize=4.3, color="#2C5E86",
                annotation_clip=False, zorder=6)
    return top


def _draw_transcripts(ax, txs, x0, x1, y0, ax_w_pt, n_max=3, dy=0.50):
    """RNA-seq 转录本轨道：外显子实心块 + 内含子细线 + 方向箭头 + 名字/TPM。"""
    sel = [t for t in txs if t["e"] >= x0 and t["s"] <= x1]
    if not sel:
        return 0, []
    sel.sort(key=lambda t: (-(t.get("tpm") or 0.0), t["s"]))
    rows = _pack(sel, gap=max(30, int((x1 - x0) * 0.01)))
    shapes = []
    for i, row in enumerate(rows[:n_max]):
        y = y0 - i * dy
        for t in row:
            col = "#3F5E86" if t["strand"] == "+" else "#8C5A3C"
            ex = t.get("exons") or [(t["s"], t["e"])]
            for a, b in ex:
                ax.add_patch(Rectangle((a, y - 0.10), max(b - a + 1, 1), 0.20,
                                       facecolor=col, edgecolor=col, lw=0.4,
                                       alpha=0.78, clip_on=True, zorder=2))
            for (a1, b1), (a2, b2) in zip(ex, ex[1:]):
                ax.plot([b1, a2], [y, y], lw=0.6, color=col, zorder=1)
            for a, b in ex:
                _chevrons(ax, a, b, y, t["strand"], ax_w_pt, step_pt=18.0,
                          size=2.4, alpha=0.85)
            _n = t.get("n_sample") or 0
            lab = "%s  %.1f%s" % (str(t["id"])[:16], t.get("tpm") or 0.0,
                                  ("  n=%d" % _n) if _n > 1 else "")
            ax.annotate(lab, xy=(max(t["s"], x0), y + 0.10), xycoords="data",
                        xytext=(1, 1), textcoords="offset points", ha="left",
                        va="bottom", fontsize=4.0, color="#33475B",
                        annotation_clip=False, zorder=6)
    return min(len(rows), n_max), shapes


def _annot_tis(ev, cds_all, strand):
    """与 ORF 重叠最多的注释 CDS → (注释起始坐标 1-based, CDS 记录)。

    正链取 CDS 起点，负链取 CDS 终点（基因 5' 端）。
    只有"注释基因校正"类事件、或 ORF 有一半以上落在该 CDS 内时才返回，
    否则基因间区/反义 ORF 会被误标一个不属于它的起始点。
    """
    s, e = int(ev["orf_start"]), int(ev["orf_end"])
    best = None
    for c in cds_all:
        if c["strand"] != strand:
            continue
        ov = min(e, c["e"]) - max(s, c["s"]) + 1
        if ov <= 0:
            continue
        if best is None or ov > best[0]:
            best = (ov, c)
    if not best:
        return None, None
    ov, c = best
    if (ev.get("event_type") != "annotated_correction"
            and ov / max(1, e - s + 1) < 0.5):
        return None, None
    return (c["s"] if strand == "+" else c["e"]), c


def _tis_codon(seq, tis, strand):
    """注释起始密码子（正链取 tis..tis+2，负链取互补）。"""
    if not seq or tis is None:
        return ""
    try:
        tri = seq[tis - 1:tis + 2] if strand == "+" else seq[tis - 3:tis]
    except Exception:
        return ""
    if len(tri) < 3:
        return ""
    tri = tri.upper()
    if strand == "-":
        tri = tri.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    return tri


def _orf_bar(ax, s, e, strand, col, tis, y, fs=5.6, codon="", ax_w_pt=600.0,
             frame_txt="", mark=None, shapes=None, ext_caption=True):
    """ORF 条：浅绿填充、黑边；名字与读框写在条正上方，条内只有箭头。"""
    ext = None
    if tis is not None:
        if strand == "+" and tis - 1 >= s:
            ext = (s, tis - 1)
        elif strand == "-" and tis + 1 <= e:
            ext = (tis + 1, e)
    body = (s, e) if not ext else ((tis, e) if strand == "+" else (s, tis))
    if ext:
        a, b = ext
        ax.add_patch(Rectangle((a, y - 0.18), b - a + 1, 0.36,
                               facecolor="#DCEFD4", edgecolor="black",
                               lw=0.6, ls=(0, (2, 1.2)), alpha=0.95,
                               clip_on=True, zorder=2))
        naa = (b - a + 1) // 3
        _x0, _x1 = ax.get_xlim()
        if ext_caption and (b - a) > (_x1 - _x0) * 0.06:
            ax.annotate("5' extension (%d aa)" % naa,
                        xy=((a + b) / 2.0, y - 0.18), xycoords="data",
                        xytext=(0, -6), textcoords="offset points",
                        ha="center", va="top", fontsize=4.6, color="#2E5E2E",
                        annotation_clip=False, zorder=6)
    ax.add_patch(Rectangle((body[0], y - 0.18), body[1] - body[0] + 1, 0.36,
                           facecolor="#CFE8C4", edgecolor="black", lw=0.6,
                           clip_on=True, zorder=2))
    if shapes is not None:
        shapes.append((body[0], body[1], y, strand))
    # 条正上方的文字（不放条内）
    lab = "novel ORF" + (("   " + frame_txt) if frame_txt else "")
    ax.annotate(lab, xy=(s + (e - s) / 2.0, y + 0.18), xycoords="data",
                xytext=(0, 2), textcoords="offset points", ha="center",
                va="bottom", fontsize=fs, color="#1F3D1F",
                annotation_clip=False, zorder=6)
    if mark is not None:
        xm = mark - (0.5 if strand == "+" else -0.5)
        ax.axvline(xm, color="#B45309", lw=0.8, ls=(0, (2.5, 1.6)),
                   clip_on=True, zorder=4)



def _enz(p):
    import re
    return set(re.findall(r"TL|LL", str(p.get("enzymes", "")).upper()))



def igv_page_v2(fig, ev, cds_all, peps, win, alt_peps=None, seq=None, clen=None,
                txs=None, models=None, depth=None):
    from sjh_pdfs import SUB_COL, EVENT_COL
    fig.clf()
    alt_peps = alt_peps or []
    ws, we = win
    s, e = int(ev["orf_start"]), int(ev["orf_end"])
    strand = ev.get("strand", "+")
    ax_w_pt = fig.get_figwidth() * 0.905 * 72.0
    tis, _hit = _annot_tis(ev, cds_all, strand)
    codon = _tis_codon(seq, tis, strand)
    # 延伸段着色只在"注释基因校正 + 与注释 CDS 同框"时成立；
    # 不同框（same_strand_other_frame）的 ORF 不是延伸，只画起始参考线。
    ext_tis = None
    if tis is not None and ev.get("event_type") == "annotated_correction":
        fo = _frame_on_contig(s, e, strand, clen)
        fc = _frame_on_contig(_hit["s"], _hit["e"], strand, clen) if _hit else None
        if fo is None or fc is None or fo == fc:
            ext_tis = tis
    # ORF 自身读框：L2 名字里自带（与建库同基准），其余按 contig 坐标算
    orf_name = str(ev.get("orf", ""))
    fr = re.search(r"_([+-])(\d)_\d+$", orf_name) if orf_name.startswith("L2_") else None
    if fr:
        orf_frame = "frame %s%s" % (fr.group(1), fr.group(2))
    else:
        orf_frame = _frame_txt(s, e, strand, clen)

    pad = max(120, int((e - s) * 0.06))
    z0, z1 = s - pad, e + pad
    if z1 - z0 > 6000:
        mid = (s + e) // 2
        z0, z1 = mid - 3000, mid + 3000

    fig.text(0.5, 0.965, "%s | %s / %s" % (ev.get("event_id", ""),
             ev.get("event_type", ""), ev.get("event_subtype", "")),
             ha="center", va="center", fontsize=10, fontweight="bold")
    fig.text(0.5, 0.930, "%s  (%s)  %s-%s" % (ev.get("orf", ""), strand,
             f"{s:,}", f"{e:,}"), ha="center", va="center", fontsize=8,
             color="#333333")
    # ---- 转录本支持（证据，不是覆盖度）
    _nm = int(ev.get("tx_models") or 0)
    _ext = str(ev.get("tx_ext") or "")
    if _nm:
        _txt = "RNA-seq: %d model%s the ORF · %s/%s samples · max TPM %s" % (
            _nm, " covers" if _nm == 1 else "s cover", ev.get("tx_samples", "?"),
            ev.get("rna_n_total", "?"), ev.get("tx_max_tpm", "?"))
        if _ext in ("yes", "no"):
            if _ext == "yes":
                _txt += " · extension covered"
                _col = "#1F7A3D"
            else:
                _dp = None
                if depth and tis:
                    _dp = (_mean_depth(depth, s, tis - 1) if strand == "+"
                           else _mean_depth(depth, tis + 1, e))
                if _dp and _dp >= 5:
                    _txt += " · extension: reads %.0f× but no transcript model" % _dp
                else:
                    _txt += " · extension NOT covered"
                _col = "#B45309"
        else:
            _col = "#1F7A3D"
    else:
        _txt = "RNA-seq: no assembled transcript overlaps this ORF"
        _col = "#6B7280"
    fig.text(0.5, 0.902, _txt, ha="center", va="bottom", fontsize=7,
             color=_col)

    col = SUB_COL.get(ev.get("event_subtype"),
                      EVENT_COL.get(ev.get("event_type"), "#888"))

    # ============ 上面板 overview
    ax1 = fig.add_axes([0.085, 0.575, 0.890, 0.320])
    _pad = (we - ws) * 0.05
    ax1.set_xlim(ws - _pad, we + _pad)
    n1, sh1 = _draw_cds(ax1, cds_all, ws, we, base_y=2.6, ax_w_pt=ax_w_pt,
                        orf_s=s - 1, orf_strand=strand, show_frame=False, clen=clen)
    y_o1 = 1.5
    so1 = []
    _orf_bar(ax1, s, e, strand, col, ext_tis, y_o1, fs=5.6, ax_w_pt=ax_w_pt,
             frame_txt=orf_frame, mark=tis, shapes=so1)
    ylo, yhi = (0.9, 2.6 + max(n1 - 1, 0) + 1.15)
    ax1.set_ylim(ylo, yhi)
    # 放大区域指示：overview 上画虚线框，标出 zoom 面板的范围
    ax1.add_patch(Rectangle((z0, ylo), z1 - z0, yhi - ylo, fill=False,
                            edgecolor="#333333", lw=0.7, ls=(0, (3, 2)),
                            clip_on=True, zorder=7))
    ax1.annotate("zoom", xy=(z0, yhi), xytext=(2, -2), textcoords="offset points",
                 ha="left", va="top", fontsize=5.0, color="#333333",
                 annotation_clip=False, zorder=7)
    _arrows(ax1, sh1 + so1, ax_w_pt, ylo, yhi,
            fig.get_figheight() * 0.320 * 72.0, step_pt=20.0)
    ax1.set_yticks([])
    _left_labels(ax1, [("overview", yhi - (yhi - ylo) * 0.06, "#33475B"),
                       ("annotated CDS", 2.6 + max(n1 - 1, 0) / 2.0, "#33475B"),
                       ("novel ORF", y_o1, "#1F3D1F")], (ylo, yhi))
    for sp in ax1.spines.values():
        sp.set_linewidth(0.7)
    ax1.tick_params(axis="x", labelsize=6.5, width=0.7, length=2.5)
    _int_ticks(ax1, ws, we)

    # ============ 下面板 zoom
    ax2 = fig.add_axes([0.085, 0.160, 0.890, 0.320])
    _pad2 = (z1 - z0) * 0.04
    ax2.set_xlim(z0 - _pad2, z1 + _pad2)
    n2, sh2 = _draw_cds(ax2, cds_all, z0, z1, base_y=4.6, ax_w_pt=ax_w_pt,
                        orf_s=s - 1, orf_strand=strand, show_frame=True, clen=clen)
    y_o2 = 3.35
    so2 = []
    _orf_bar(ax2, s, e, strand, col, ext_tis, y_o2, fs=5.8, ax_w_pt=ax_w_pt,
             codon=codon, frame_txt=orf_frame, mark=tis, shapes=so2,
             ext_caption=not (txs or []))

    # ---- 肽段两行 + 镜像核心阴影（同座位别名 ORF 的肽用淡色，不编号）
    own_seqs = {str(p.get("peptide", "")) for p in peps}
    alt_peps = [p for p in alt_peps
                if str(p.get("peptide", "")) not in own_seqs]
    t_peps = [p for p in peps if "TL" in _enz(p)]
    l_peps = [p for p in peps if "LL" in _enz(p)]
    t_alt = [p for p in alt_peps if "TL" in _enz(p)]
    l_alt = [p for p in alt_peps if "LL" in _enz(p)]
    y_t, y_l = 0.50, 0.10
    def _pep_rows(plist, y, color, alpha=0.95):
        for p in sorted(plist, key=lambda q: -(int(q["pep_nt_end"])
                                               - int(q["pep_nt_start"]))):
            try:
                ps, pe = int(p["pep_nt_start"]), int(p["pep_nt_end"])
            except (KeyError, TypeError, ValueError):
                continue
            if pe < z0 or ps > z1:
                continue
            ax2.add_patch(Rectangle((ps, y - 0.11), max(pe - ps + 1, 1), 0.22,
                                    facecolor=color, edgecolor="black", lw=0.5,
                                    alpha=alpha, clip_on=True, zorder=3))
        return len([p for p in plist if z0 <= int(p["pep_nt_end"])
                    and int(p["pep_nt_start"]) <= z1])

    _pep_rows(t_alt, y_t, "#2F6FBA", alpha=0.5)
    _pep_rows(l_alt, y_l, "#C86A3F", alpha=0.5)
    _pep_rows(t_peps, y_t, "#2F6FBA")
    _pep_rows(l_peps, y_l, "#C86A3F")

    n_d = _draw_depth(ax2, depth or [], z0, z1, y_base=2.55) if depth else 0.0
    n_tx, _ = _draw_transcripts(ax2, txs or [], z0, z1, y0=2.05, ax_w_pt=ax_w_pt)
    ylo2, yhi2 = (y_l - 0.42, 4.6 + max(n2 - 1, 0) + 1.15)
    ax2.set_ylim(ylo2, yhi2)
    _arrows(ax2, sh2 + so2, ax_w_pt, ylo2, yhi2,
            fig.get_figheight() * 0.320 * 72.0, step_pt=15.0)
    ax2.set_yticks([])
    _lab2 = [("zoom", yhi2 - (yhi2 - ylo2) * 0.06, "#33475B"),
             ("annotated CDS", 4.6 + max(n2 - 1, 0) / 2.0, "#33475B"),
             ("novel ORF", y_o2, "#1F3D1F"),
             ("Trypsin", y_t, "#2F6FBA"),
             ("LysargiNase", y_l, "#C86A3F")]
    if n_d:
        _lab2.insert(3, ("RNA-seq" + chr(10) + "depth", 2.55 + 0.22, "#2C5E86"))
    if n_tx:
        _lab2.insert(4 if n_d else 3,
                     ("StringTie" + chr(10) + "transcripts",
                      2.05 - (n_tx - 1) * 0.21, "#3F5E86"))
    _left_labels(ax2, _lab2, (ylo2, yhi2))

    # ---- 注释列（右侧竖虚线 + 序号，T 先 L 后、组内 N→C）
    idx = 0
    for plist, y in ((t_peps, y_t), (l_peps, y_l)):
        for p in sorted([q for q in plist
                         if int(q.get("pep_nt_end") or 0) >= z0
                         and int(q.get("pep_nt_start") or 0) <= z1],
                        key=lambda q: int(q.get("pep_nt_start") or 0)):
            idx += 1
            ps, pe = int(p["pep_nt_start"]), int(p["pep_nt_end"])
            ax2.text((ps + pe) / 2.0, y + 0.13, str(idx), ha="center",
                     va="bottom", fontsize=4.8, color="#374151", zorder=6)

    ax2.set_xlabel("%s position (bp)" % (ev.get("seqid") or "genome"),
                   fontsize=7.5)
    for sp in ax2.spines.values():
        sp.set_linewidth(0.7)
    ax2.tick_params(axis="x", labelsize=6.5, width=0.7, length=2.5)
    _int_ticks(ax2, z0, z1)

    # ---- RNA-seq 丰度（有值才画第二行）

    fig.text(0.5, 0.092, "page %s | n_peptide %s | n_psm %s | tier %s | "
             "annotated CDS %s | TIS %s (dashed line) | alias peptides %d (pale, not in n_peptide) | window %s-%s"
             % (ev.get("_page", "?"), ev.get("n_peptide", ""),
                ev.get("n_psm", ""), ev.get("tier_mix", ""),
                len([c for c in cds_all if c["e"] >= z0 and c["s"] <= z1]),
                (("%s %s" % (f"{tis:,}", codon)).strip() if tis else "n.a."),
                len(alt_peps),
                f"{z0:,}", f"{z1:,}"), ha="center", va="center", fontsize=6,
             color="#6B7280")
