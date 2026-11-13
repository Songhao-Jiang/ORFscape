"""sjh_pdfs.py -- rebuild the two SJH deliverables in the isoPSV native layout.

1) isoPSV_spectra.pdf  : A4 page per novel peptide, laid out exactly like
   isoPSV's `novel_peptide_spectra.pdf` (header block, ion ladder with b/y
   subscripts and cleavage ticks, grey experimental sticks + coloured
   theoretical ions, Base peak / Mass error, prediction band, ppm band,
   two candidate panels).
2) IGV_like_annotation.pdf : one page per gene event, tracks packed inside a
   framed axes, overlapping CDS on separate rows, peptides as rectangles.

All file paths are CLI parameters; nothing is joined or guessed here.
"""
import argparse
import csv
import math
import os
import pickle
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.markers import CARETLEFT, CARETRIGHT
from matplotlib.backends.backend_pdf import PdfPages

# ---------------------------------------------------------------- constants
AA_MASS = {
    "A": 71.0371103, "C": 103.0091803, "D": 115.0269385, "E": 129.0425877,
    "F": 147.0684087, "G": 57.0214611, "H": 137.0589059, "I": 113.0840579,
    "K": 128.0949557, "L": 113.0840579, "M": 131.0404787, "N": 114.0429222,
    "P": 97.0527595, "Q": 128.0585714, "R": 156.1011021, "S": 87.0320244,
    "T": 101.0476736, "V": 99.0684087, "W": 186.0793065, "Y": 163.0633228,
    "X": 113.0840579, "U": 150.9536356,
}
MOD_MASS = {
    "Oxidation[M]": 15.994915, "Carbamidomethyl[C]": 57.021464,
    "Acetyl[AnyN-term]": 42.010565, "Acetyl[ProteinN-term]": 42.010565,
    "Acetyl[K]": 42.010565, "Deamidated[N]": 0.984016, "Deamidated[Q]": 0.984016,
    "Gln->pyro-Glu[AnyN-termQ]": -17.026549, "Glu->pyro-Glu[AnyN-termE]": -18.010565,
    "Amidated[AnyC-term]": -0.984016, "Methyl[AnyC-term]": 14.01565,
    "Carbamyl[AnyN-term]": 43.005814, "Formyl[AnyN-term]": 27.994915,
    "Dimethyl[K]": 28.0313, "Trimethyl[K]": 42.04695,
    "Phospho[S]": 79.966331, "Phospho[T]": 79.966331, "Phospho[Y]": 79.966331,
    "Carboxymethyl[C]": 58.005479,
}
PROTON = 1.007276
H2O = 18.010565
CO = 27.994915
NH3 = 17.026549

ION_COL = {"a": "#004000", "b": "#008000", "c": "#CA1D52", "x": "#650E29",
           "y": "#CA1D52", "z": "#943AA4", "M": "#808080",
           "I": "#FF8C00", "O": "#9ACD32"}
ION_PRI = {"b": 1, "y": 2, "a": 3, "c": 4, "x": 5, "z": 6, "M": 7, "I": 8, "O": 9}
TOL_PPM = 20.0


def iln(s):
    return (s or "").upper().replace("I", "L")


def short_orf(orf):
    """压缩超长 ORF 显示标签（pAnno 人基因组版带过来的修正）。

    转录本 ORF 的前缀可能极长（rnaSPAdes 全路径 / GCF 长 accession），
    图上放不下；取 SRR 号 + 末尾基因段 + 坐标，短名原样返回。
    """
    s = str(orf or "")
    if ":" not in s:
        return s if len(s) <= 60 else s[:57] + "..."
    prefix, spec = s.rsplit(":", 1)
    if prefix.startswith("GCF_") and ".fasta_line" in prefix:
        seqid = prefix.split(".fasta_line", 1)[1].split("_", 1)[1]
    elif prefix.startswith("GCF_"):
        seqid = prefix.rsplit("_", 1)[-1]
    else:
        parts = prefix.split("_")
        gene = "_".join(parts[-2:]) if len(parts) >= 2 else prefix
        m = re.search(r"SRR\d+", prefix)
        seqid = f"{m.group(0)}_{gene}" if m else gene
    return "%s:%s" % (seqid, spec)


def split_mods(mod):
    """'21,Oxidation[M];22,Oxidation[M];' -> [(21,'Oxidation[M]'), ...]"""
    out = []
    for tok in (mod or "").split(";"):
        tok = tok.strip()
        if not tok:
            continue
        p, _, name = tok.partition(",")
        try:
            out.append((int(p), name.strip()))
        except ValueError:
            continue
    return out


def peptide_mass(seq, mod):
    """Return (residue neutral masses with mods applied, modified positions)."""
    res = [AA_MASS.get(a, 0.0) for a in seq]
    pos = []
    for p, name in split_mods(mod):
        dm = MOD_MASS.get(name)
        if dm is None:
            dm = MOD_MASS.get(re.sub(r"\[[^\]]*\]$", "", name) + "]")
        dm = 0.0 if dm is None else dm
        if p <= 0:
            res[0] += dm
            pos.append(1)
        elif p <= len(res):
            res[p - 1] += dm
            pos.append(p)
    return res, sorted(set(pos))


def theoretical_ions(seq, mod):
    """b/a/c/y/x/z (1+/2+) plus the M series, as (mz, label, ion, num, chg)."""
    if not seq:
        return []
    res, _ = peptide_mass(seq, mod)
    n = len(res)
    cum = np.cumsum(res)
    total = cum[-1]
    out = []
    for k in range(1, n):
        neu_b = cum[k - 1]
        neu_y = total - cum[n - k - 1]
        out.append((neu_b + PROTON, "b%d+" % k, "b", k, 1))
        out.append((neu_b - CO + PROTON, "a%d+" % k, "a", k, 1))
        out.append((neu_b + NH3 + PROTON, "c%d+" % k, "c", k, 1))
        out.append(((neu_b + 2 * PROTON) / 2.0, "b%d++" % k, "b", k, 2))
        out.append(((neu_b - CO + 2 * PROTON) / 2.0, "a%d++" % k, "a", k, 2))
        out.append(((neu_b + NH3 + 2 * PROTON) / 2.0, "c%d++" % k, "c", k, 2))
        y1 = neu_y + H2O
        out.append((y1 + PROTON, "y%d+" % k, "y", k, 1))
        out.append((y1 + CO + PROTON, "x%d+" % k, "x", k, 1))
        out.append((y1 - NH3 + PROTON, "z%d+" % k, "z", k, 1))
        out.append(((y1 + 2 * PROTON) / 2.0, "y%d++" % k, "y", k, 2))
        out.append(((y1 + CO + 2 * PROTON) / 2.0, "x%d++" % k, "x", k, 2))
        out.append(((y1 - NH3 + 2 * PROTON) / 2.0, "z%d++" % k, "z", k, 2))
    pm = total + H2O
    mh = pm + PROTON
    out.append((mh, "M+", "M", 0, 1))
    out.append((mh - H2O, "M-H2O+", "M", 0, 1))
    out.append((mh - NH3, "M-NH3+", "M", 0, 1))
    out.append(((mh + PROTON) / 2.0, "M++", "M", 0, 2))
    out.append(((mh + PROTON - H2O) / 2.0, "M-H2O++", "M", 0, 2))
    out.append(((mh + PROTON - NH3) / 2.0, "M-NH3++", "M", 0, 2))
    return out, mh


def match_ions(ions, mz, ints, tol_ppm=TOL_PPM):
    """Attach the strongest peak inside +/- tol to every theoretical ion."""
    if len(mz) == 0:
        return [(m, lab, ion, num, chg, None, 0.0) for (m, lab, ion, num, chg) in ions]
    order = np.argsort(mz)
    mz_s = mz[order]
    ints_s = ints[order]
    res = []
    for (m, lab, ion, num, chg) in ions:
        lo = m * (1.0 - tol_ppm * 1e-6)
        hi = m * (1.0 + tol_ppm * 1e-6)
        i0 = np.searchsorted(mz_s, lo, side="left")
        i1 = np.searchsorted(mz_s, hi, side="right")
        if i1 <= i0:
            res.append((m, lab, ion, num, chg, None, 0.0))
            continue
        seg = ints_s[i0:i1]
        j = i0 + int(np.argmax(seg))
        res.append((m, lab, ion, num, chg, float(mz_s[j]), float(ints_s[j])))
    return res


# ------------------------------------------------------------------ page io
def load_table(path):
    with open(path, encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def best_psms(spectra_paths, wanted):
    """Map iln(sequence) -> the highest Final_Score PSM row."""
    best = {}
    for p in spectra_paths:
        with open(p, encoding="utf-8-sig", errors="replace", newline="") as fh:
            rd = csv.DictReader(fh, delimiter="\t")
            for row in rd:
                key = iln(row.get("Sequence", ""))
                if key not in wanted:
                    continue
                try:
                    sc = float(row.get("Final_Score") or 0.0)
                except ValueError:
                    sc = 0.0
                cur = best.get(key)
                if cur is None or sc > cur[0]:
                    best[key] = (sc, row)
    return {k: v[1] for k, v in best.items()}


def split_title(fn):
    """XP_..._T5_45.15665.15665.4.0.dta -> (stem, '15665')"""
    base = fn.split("/")[-1].split("\\")[-1]
    if base.lower().endswith(".dta"):
        base = base[:-4]
    parts = base.split(".")
    if len(parts) >= 4:
        return parts[0], parts[1]
    return base, ""


def read_predictions(path):
    """pDeep2/ms2pip predict.txt -> {(peptide, charge): {series: [intensities]}}"""
    out = {}
    seq = None
    chg = 0
    series = {}
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("BEGIN IONS"):
                seq, chg, series = None, 0, {}
            elif line.startswith("TITLE="):
                t = line.split("=", 1)[1]
                parts = t.split("|")
                seq = parts[0]
                for q in parts[1:]:
                    if q.isdigit():
                        chg = int(q)
            elif "=" in line:
                k, _, v = line.partition("=")
                if k in ("b+1", "b+2", "y+1", "y+2"):
                    try:
                        series[k] = [float(x) for x in v.split(",") if x]
                    except ValueError:
                        pass
            elif line.startswith("END IONS"):
                if seq:
                    out[(seq, chg)] = series
    return out


def theory_cosine(mz, ints, sticks, tol_ppm=TOL_PPM):
    """Cosine between the measured spectrum and a predicted ion series.

    Each predicted b/y ion is paired with the strongest measured peak inside
    the +/- tol window (no peak -> 0); both vectors are sqrt-scaled, which is
    the usual scale for predicted-intensity comparison.
    """
    if not sticks or len(mz) == 0:
        return None
    order = np.argsort(mz)
    mzs, ins = mz[order], ints[order]
    a, b = [], []
    for (m, v, _ion) in sticks:
        lo, hi = m * (1 - tol_ppm * 1e-6), m * (1 + tol_ppm * 1e-6)
        i0 = np.searchsorted(mzs, lo, side="left")
        i1 = np.searchsorted(mzs, hi, side="right")
        a.append(v)
        b.append(float(ins[i0:i1].max()) if i1 > i0 else 0.0)
    a = np.sqrt(np.asarray(a, float))
    b = np.asarray(b, float)
    if b.max() <= 0 or a.max() <= 0:
        return 0.0
    b = np.sqrt(b / b.max())
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    return float(a @ b / (na * nb)) if na > 0 and nb > 0 else 0.0


def pred_sticks(seq, mod, series):
    """Predicted b/y intensities -> [(mz, height 0..1, ion)]"""
    if not series or not seq:
        return []
    res, _ = peptide_mass(seq, mod)
    n = len(res)
    cum = np.cumsum(res)
    total = cum[-1]
    sticks = []
    for k in range(1, n):
        for key, ion in (("b+1", "b"), ("b+2", "b"), ("y+1", "y"), ("y+2", "y")):
            arr = series.get(key)
            if not arr or len(arr) < n - 1:
                continue
            v = arr[k - 1]
            if v <= 0:
                continue
            if ion == "b":
                neu = cum[k - 1]
                mz = neu + PROTON if key == "b+1" else (neu + 2 * PROTON) / 2.0
            else:
                neu = total - cum[n - k - 1] + H2O
                mz = neu + PROTON if key == "y+1" else (neu + 2 * PROTON) / 2.0
            sticks.append((mz, v, ion))
    return sticks


# --------------------------------------------------------------- isoPSV page
def draw_ladder(ax, seq, modpos, ions, charge):
    """Sequence row above the spectrum: letters, b/y subscripts, cleavage ticks."""
    if not seq:
        return
    n = len(seq)
    pos = ax.get_position()
    axw_in = pos.width * ax.figure.get_figwidth()
    axw_in = axw_in if axw_in > 0.5 else 7.0
    axh_in = pos.height * ax.figure.get_figheight()
    axh_in = axh_in if axh_in > 0.2 else 4.0
    x0, x1 = 0.06, 0.94                      # axes-fraction span of the ladder
    slot = min(0.16 / axw_in, (x1 - x0) / n)   # original isoPSV ladder pitch
    xc = (x0 + x1) / 2.0 - n * slot / 2.0
    xs = [xc + (i + 0.5) * slot for i in range(n)]
    ymid = 1.0 + 0.260 / axh_in          # all offsets in inches, not axes fractions
    tick = 0.075 / axh_in                # must clear the 8-pt glyphs on both sides
    b_sites, y_sites = [], []
    for (m, lab, ion, num, chg, pk, inten) in ions:
        if pk is None:
            continue
        if ion == "b":
            b_sites.append(num)
        elif ion == "y":
            y_sites.append(n - num)
    b_sites = sorted(set(b_sites))
    y_sites = sorted(set(y_sites))
    sites = sorted(set(b_sites) | set(y_sites))
    for s in sites:
        if s >= n:
            continue
        vx = xs[s - 1] + slot * 0.5          # cleavage sits BETWEEN residue s and s+1
        ax.plot([vx, vx], [ymid - tick, ymid + tick], color="black", lw=0.5,
                transform=ax.transAxes, clip_on=False, solid_capstyle="butt")
    slot_pt = slot * axw_in * 72.0
    fs = max(5.0, min(8.0, slot_pt * 0.70))
    fl = max(4.0, min(6.0, slot_pt * 0.52))
    thin = slot_pt < 10.0                    # only very tight ladders thin labels
    for k, s in enumerate(b_sites):
        if s < n and not (thin and k % 2):
            vx = xs[s - 1] + slot * 0.5
            ax.plot([vx, vx - slot * 0.5], [ymid - tick, ymid - tick],
                    color="black", lw=0.5, transform=ax.transAxes,
                    clip_on=False, solid_capstyle="butt")
            ax.text(xs[s - 1], ymid - 0.125 / axh_in, r"$b_{%d}$" % s,
                    transform=ax.transAxes, ha="center", va="center",
                    fontsize=fl, color=ION_COL["b"], clip_on=False)
    for k, s in enumerate(y_sites):
        if s < n and not (thin and k % 2):
            vx = xs[s - 1] + slot * 0.5
            ax.plot([vx, vx + slot * 0.5], [ymid + tick, ymid + tick],
                    color="black", lw=0.5, transform=ax.transAxes,
                    clip_on=False, solid_capstyle="butt")
            ax.text(xs[s], ymid + 0.125 / axh_in, r"$y_{%d}$" % (n - s),
                    transform=ax.transAxes, ha="center", va="center",
                    fontsize=fl, color=ION_COL["y"], clip_on=False)
    for i, ch in enumerate(seq):
        res_i = i + 1
        bold = "bold" if res_i in modpos else "normal"
        col = ION_COL["I"] if res_i in modpos else "black"
        ax.text(xs[i], ymid, ch, transform=ax.transAxes, ha="center",
                va="center", fontsize=fs, fontweight=bold, color=col,
                clip_on=False)
    if charge:
        ax.text(xs[0] - slot * 1.6, ymid, "%d+" % charge, transform=ax.transAxes,
                ha="right", va="center", fontsize=min(fs, 6.8), fontweight="bold",
                color="#CA1D52", clip_on=False)


def panel_spectrum(ax, mz, ints_n, ions, head, ann, charge, seq, modpos, xmax):
    ax.set_xlim(0, xmax)
    ax.set_ylim(0, 115)
    ax.set_yticks([0, 20, 40, 60, 80, 100])
    ax.set_yticklabels(["0", "20", "40", "60", "80", "100"], fontsize=8)
    ax.set_xticks([])
    ax.set_ylabel("Relative intensity (%)", fontsize=8.5)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.vlines(mz, 0, ints_n * 100, lw=0.35, color="#999999")

    hit = [t for t in ions if t[5] is not None]
    if hit:
        hit.sort(key=lambda t: (ION_PRI.get(t[2], 9), t[3], t[4]))
        groups = {}
        for t in hit:
            groups.setdefault(t[5], []).append(t)
        for pk, lst in groups.items():
            lst.sort(key=lambda t: (ION_PRI.get(t[2], 9), t[3], t[4]))
            h = float(np.interp(pk, mz, ints_n * 100)) if len(mz) else 0.0
            ax.vlines(pk, 0, h, lw=0.9, color=ION_COL[lst[0][2]])
            ax.text(pk, h + 2.0, "".join(t[1] for t in lst), rotation=90,
                    fontsize=6.0, color=ION_COL[lst[0][2]], ha="center",
                    va="bottom", clip_on=True)

    for i, line in enumerate(ann):
        ax.text(0.995, 0.935 - 0.070 * i, line, transform=ax.transAxes,
                ha="right", va="top", fontsize=8.0, clip_on=False)
    if head:
        hf = 7.0 if len(head) <= 95 else (6.4 if len(head) <= 115 else 5.8)
        ahi = ax.get_position().height * ax.figure.get_figheight()
        ax.text(0.0, 1.0 + 0.510 / (ahi if ahi > 0.2 else 4.0), head,
                transform=ax.transAxes, ha="left", va="bottom",
                fontsize=hf, clip_on=False)
    draw_ladder(ax, seq, modpos, ions, charge)


def panel_pred(ax, sticks, xmax):
    ax.set_xlim(0, xmax)
    ax.set_ylim(0, 100)
    ax.set_xticks([])
    ax.set_yticks([20, 60, 100])
    ax.set_yticklabels(["20", "60", "100"], fontsize=6.2)
    ax.tick_params(axis="y", length=2, width=0.6, pad=1.5)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    if sticks:
        for (m, v, ion) in sticks:
            ax.vlines(m, 0, v * 100, lw=0.7, color=ION_COL.get(ion, "#999999"))
        ax.text(0.985, 0.10, "Prediction", transform=ax.transAxes, ha="right",
                va="bottom", fontsize=7.5, color="black")
    else:
        ax.text(0.5, 0.5, "no prediction", transform=ax.transAxes, ha="center",
                va="center", fontsize=9.5, color="#777777")
        ax.text(0.985, 0.10, "Prediction", transform=ax.transAxes, ha="right",
                va="bottom", fontsize=7.5, color="black")


def panel_ppm(ax, ions, mz, xmax, show_xlabel=True):
    ax.set_xlim(0, xmax)
    ax.set_ylim(-20, 20)
    ax.set_yticks([-20, 0, 20])
    ax.set_yticklabels(["-20", "0", "20"], fontsize=8)
    ax.set_ylabel("ppm", fontsize=8.5)
    ax.axhline(0, color="#999999", lw=0.6)
    ax.axhline(20, color="#BBBBBB", lw=0.6, ls=(0, (3, 3)))
    ax.axhline(-20, color="#BBBBBB", lw=0.6, ls=(0, (3, 3)))
    for (m, lab, ion, num, chg, pk, inten) in ions:
        if pk is None:
            continue
        err = (pk - m) / m * 1e6
        if -20 <= err <= 20:
            ax.plot(pk, err, "o", ms=2.6, color=ION_COL.get(ion, "#999999"),
                    mew=0)
    ax.tick_params(axis="x", labelsize=8)
    if show_xlabel:
        ax.set_xlabel("m/z", fontsize=8.5, style="italic")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def panel_theory(ax, sticks, xmax, cos_txt, show_x=True):
    """AlphaPeptDeep spectrum, y axis inverted, same height as the experiment."""
    ax.set_xlim(0, xmax)
    ax.set_ylim(100, 0)                       # inverted: intensity grows downward
    ax.set_yticks([20, 60, 100])
    ax.set_yticklabels(["20", "60", "100"], fontsize=8)
    ax.set_ylabel("Predicted (%)", fontsize=8.0, labelpad=2)
    ax.tick_params(axis="x", labelsize=8)
    if show_x:
        ax.set_xlabel("m/z", fontsize=8.5, style="italic")
    else:
        ax.tick_params(axis="x", labelbottom=False)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    vis = [t for t in sticks if 0 <= t[0] <= xmax]
    if vis:
        vmax = max(v for (_m, v, _i) in vis) or 1.0   # scale to the visible window
        for (m, v, ion) in vis:
            ax.vlines(m, 0, v / vmax * 100.0, lw=0.9,
                      color=ION_COL.get(ion, "#999999"))
    else:
        ax.text(0.5, 0.5, "no prediction", transform=ax.transAxes, ha="center",
                va="center", fontsize=9.5, color="#777777")


def iso_page(fig, row, psm, peaks, pred, xmax_hint=None):
    fig.clf()
    seq = row["peptide"]
    mod = row["modification"]
    modpos = sorted({p for p, _ in split_mods(mod)})
    mz = np.array([p[0] for p in peaks], dtype=float) if peaks else np.array([])
    it = np.array([p[1] for p in peaks], dtype=float) if peaks else np.array([])
    if len(mz):
        o = np.argsort(mz)
        mz, it = mz[o], it[o]
        base = np.max(it) if np.max(it) > 0 else 1.0
    else:
        base = 1.0
    ints_n = it / base if len(it) else it
    xmax = (np.max(mz) * 1.03) if len(mz) else 1000.0

    ions_t, calc_mh = theoretical_ions(seq, mod)
    ions = match_ions(ions_t, mz, it)
    try:
        exp_mh = float(psm.get("Exp.MH+") or 0.0)
    except ValueError:
        exp_mh = 0.0
    try:
        charge = int(float(psm.get("Charge") or 0))
    except ValueError:
        charge = 0
    bp_lab = "%.2E" % base
    bp_lab = bp_lab.replace("E+0", "E").replace("E-0", "E-")
    ann = ["Experiment  Base peak = %s" % bp_lab]
    if exp_mh > 0 and calc_mh > 0:
        d = exp_mh - calc_mh
        ann.append("Mass error = %+.4f Da / %+.2f ppm" % (d, d / calc_mh * 1e6))
    else:
        ann.append("Mass error = -")

    prot = (psm.get("Proteins") or "").strip()
    acc = short_orf(prot.split("/")[0]) if prot and prot != "-" else "-"
    pos = (psm.get("Positions") or "").strip()
    pre, post = "-", "-"
    if pos and "[" in pos and "]" in pos:
        pre = pos.split("[")[0] or "-"
        post = pos.split("]")[1] or "-"
    elif "," in pos:
        f = [x for x in pos.replace("/", "").split(",") if x]
        if len(f) >= 2:
            pre, post = f[1] or "-", f[-1] or "-"
    head = "| %s | (%s)%s(%s) | %s |" % (acc, pre, seq, post, mod or "-")

    handle = psm.get("_basename") or psm.get("File_Name", "")
    modtxt = mod if mod else "-"
    sticks = pred_sticks(seq, mod, pred) if pred else []
    cos = theory_cosine(mz, it, sticks)
    cos_txt = "%.3f" % cos if cos is not None else "-"

    lines = [
        "Spectrum: %s" % handle,
        "Novel peptide: %s | Modification: %s" % (seq, modtxt),
        "Judgment: %s" % row.get("_judgment", "-"),
        "AlphaPeptDeep theoretical similarity: %s" % cos_txt,
    ]
    for i, ln in enumerate(lines):
        fig.text(0.5, 0.955 - i * 0.025, ln, ha="center", va="center",
                 fontsize=7.4)

    ax_sp = fig.add_axes([0.085, 0.500, 0.855, 0.250])
    ax_th = fig.add_axes([0.085, 0.250, 0.855, 0.250])
    ax_pp = fig.add_axes([0.085, 0.115, 0.855, 0.105])
    panel_spectrum(ax_sp, mz, ints_n, ions, head, ann, charge, seq, modpos, xmax)
    panel_theory(ax_th, sticks, xmax, cos_txt, show_x=False)
    panel_ppm(ax_pp, ions, mz, xmax)


# ------------------------------------------------------------------ IGV page
def parse_gtf(path):
    cds = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[2] != "CDS":
                continue
            at = {}
            for kv in p[8].split(";"):
                kv = kv.strip()
                if not kv:
                    continue
                k, _, v = kv.partition(" ")
                at[k] = v.strip().strip('"')
            locus = at.get("locus_tag") or at.get("protein_id") or ""
            # 保留 seqid：多 contig 时不能只靠数值窗口匹配注释 CDS
            cds.append(dict(seqid=p[0], s=int(p[3]), e=int(p[4]), strand=p[6],
                            label=locus, locus=locus, gene=at.get("gene", ""),
                            product=at.get("product", "")))
    return cds


def pack_rows(items, gap):
    """Greedy interval packing -> list of rows, each a list of items."""
    rows = []
    for it in sorted(items, key=lambda x: x["s"]):
        for row in rows:
            if it["s"] > row[-1]["e"] + gap:
                row.append(it)
                break
        else:
            rows.append([it])
    return rows


JUDGMENT = {
    "Tier1": "Mirror peptides",
    "Tier2": "ORF with both protease evidence",
    "Tier3": "only one protease evidence",
}

EVENT_COL = {"annotated_correction": "#C0504D",
             "intergenic_new_gene": "#4472C4"}
SUB_COL = {"nterm_extension": "#C0504D", "intergenic_only": "#4472C4",
           "same_strand_other_frame": "#E8A33D", "unmapped": "#999999"}


def igv_page(fig, ev, cds_all, peps, win):
    fig.clf()
    ws, we = win
    pad = max(50, int((we - ws) * 0.015))
    x0, x1 = ws - pad, we + pad

    # 只取事件所在 contig 的注释基因；事件没给 seqid（单 contig 项目）时回落为
    # "全部 CDS 若只剩一个 contig 就用它"，避免多 contig 数据靠数值窗口误匹配。
    _sq = ev.get("seqid") or ""
    if not _sq:
        _ids = {c.get("seqid", "") for c in cds_all if c.get("seqid")}
        if len(_ids) == 1:
            _sq = next(iter(_ids))
    sel = [c for c in cds_all
           if c["e"] >= x0 and c["s"] <= x1
           and (not _sq or c.get("seqid", "") == _sq)]
    rows = pack_rows(sel, gap=max(20, int((x1 - x0) * 0.004))) if sel else []
    n_row = max(1, len(rows))

    fig.text(0.5, 0.945, "%s | %s / %s" % (short_orf(ev["event_id"]), ev["event_type"],
             ev["event_subtype"]), ha="center", va="center", fontsize=10,
             fontweight="bold")
    fig.text(0.5, 0.905, "%s  (%s)  %s-%s" % (short_orf(ev["orf"]), ev["strand"],
             f"{int(ev['orf_start']):,}", f"{int(ev['orf_end']):,}"),
             ha="center", va="center", fontsize=8, color="#333333")

    ax = fig.add_axes([0.115, 0.115, 0.858, 0.755])
    ax.set_xlim(x0, x1)
    y_orf = -0.75
    y_pep = -1.75
    ax.set_ylim(-2.20, n_row + 0.28)
    ax.set_facecolor("#FBFBFB")
    for s in ax.spines.values():
        s.set_visible(True)
        s.set_color("#888888")
        s.set_linewidth(0.8)

    # annotated CDS, one row per non-overlapping group
    span = x1 - x0
    ax_w_pt = fig.get_figwidth() * 0.858 * 72.0     # axes width in points

    def est_w(txt, fspt):
        """Rough text width in data units (0.62 em per char, italic-safe)."""
        return len(txt) * fspt * 0.62 / ax_w_pt * span

    for i, row in enumerate(rows):
        y = n_row - 1 - i
        row.sort(key=lambda c: c["s"])
        placed = []                                  # (level, x0c, x1c) per label
        for c in row:
            w = c["e"] - c["s"] + 1
            ax.add_patch(Rectangle((c["s"], y - 0.19), w, 0.38,
                                   facecolor="#8FA3BF", edgecolor="#5A6B85",
                                   lw=0.5, clip_on=True, zorder=2))
            # 方向箭头：CDS 框内画 chevron（正链向右、负链向左），宽度按点数换算
            step = max(16.0 / ax_w_pt * span, span * 1e-5)
            n_ch = int(w / step)
            if n_ch >= 1:
                xs = [c["s"] + step * (k + 0.5) for k in range(n_ch)]
                ax.plot(xs, [y] * len(xs), linestyle="none",
                        marker=CARETRIGHT if c["strand"] == "+" else CARETLEFT,
                        markersize=2.8, markeredgewidth=0.5,
                        color="white", alpha=0.9, clip_on=True, zorder=3)
            l1 = (c["locus"] + ("  " + c["gene"] if c["gene"] else "")).strip()
            prod = c["product"]
            if len(prod) > 38:
                prod = prod[:36] + ".."
            l2 = prod
            lw = max(est_w(l1, 5.8), est_w(l2, 4.8))
            pad = span * 0.008
            cx = min(max(c["s"] + w / 2.0, x0 + lw / 2 + pad), x1 - lw / 2 - pad)
            lv = 0
            for k in (0, 1, 2):                      # stagger when labels would collide
                clash = any(p[0] == k and not (cx + lw / 2 < p[1] or cx - lw / 2 > p[2])
                            for p in placed)
                if not clash:
                    lv = k
                    break
            placed.append((lv, cx - lw / 2, cx + lw / 2))
            yb = y + 0.34 + lv * 0.30
            ax.text(cx, yb, l1, rotation=0, fontsize=5.8, ha="center", va="bottom",
                    color="#1F2A38", clip_on=True, zorder=3)
            if l2:
                ax.text(cx, yb - 0.135, l2, rotation=0, fontsize=4.8, ha="center",
                        va="bottom", color="#5A6B85", style="italic",
                        clip_on=True, zorder=3)

    ax.text(x0 - (x1 - x0) * 0.012, (n_row - 1) / 2.0,
            "annotated CDS", ha="right", va="center", fontsize=7.5,
            color="#33475B")

    # novel ORF
    s, e = int(ev["orf_start"]), int(ev["orf_end"])
    col = SUB_COL.get(ev["event_subtype"], EVENT_COL.get(ev["event_type"], "#888"))
    ax.add_patch(Rectangle((s, y_orf - 0.17), e - s + 1, 0.34, facecolor=col,
                           edgecolor=col, alpha=0.85, clip_on=True, zorder=2))
    # novel ORF 同样标方向（框内 chevron）
    _st = ">" if ev.get("strand") == "+" else "<"
    _step = max(14.0 / ax_w_pt * (x1 - x0), (x1 - x0) * 1e-5)
    _n = int((e - s) / _step)
    if _n >= 1:
        _xs = [s + _step * (k + 0.5) for k in range(_n)]
        ax.plot(_xs, [y_orf] * len(_xs), linestyle="none",
                marker=CARETRIGHT if _st == ">" else CARETLEFT,
                markersize=2.6, color="white", alpha=0.9, clip_on=True, zorder=3)
    ax.text(s + (e - s) / 2, y_orf + 0.32, short_orf(ev["orf"]),
            ha="center", va="bottom", fontsize=6.4, color=col, clip_on=True)
    ax.text(x0 - (x1 - x0) * 0.012, y_orf, "novel ORF", ha="right", va="center",
            fontsize=7.5, color="#33475B")

    # peptides as rectangles
    for p in peps:
        ps, pe = int(p["pep_nt_start"]), int(p["pep_nt_end"])
        if pe < x0 or ps > x1:
            continue
        ax.add_patch(Rectangle((ps, y_pep - 0.18), max(pe - ps + 1, 1), 0.36,
                               facecolor="#2F6FBA", edgecolor="#1B4C86",
                               lw=0.4, clip_on=True, zorder=2))
    ax.text(x0 - (x1 - x0) * 0.012, y_pep, "peptides", ha="right", va="center",
            fontsize=7.5, color="#33475B")

    # axis
    ax.set_yticks([])
    span = x1 - x0
    step = 10 ** int(math.floor(math.log10(max(span / 6.0, 1))))
    for mult in (1, 2, 5, 10):
        if span / (step * mult) <= 7:
            step *= mult
            break
    ticks = [t for t in range(int(x0 // step) * step, int(x1) + step, step)
             if x0 <= t <= x1]
    ax.set_xticks(ticks)
    ax.set_xticklabels(["{:,}".format(t) for t in ticks], fontsize=7.5)
    ax.tick_params(axis="x", length=3, width=0.7, colors="#333333")
    ax.set_xlabel(("%s position (bp)" % _sq) if _sq else "transcript position (nt)",
                  fontsize=8)
    for t in ticks[1:]:
        ax.axvline(t, color="#E4E4E4", lw=0.6, zorder=0)

    fig.text(0.5, 0.035, "page %s | n_peptide %s | n_psm %s | %s | %s | window %s-%s"
             % (ev.get("_page", "?"), ev.get("n_peptide", ""), ev.get("n_psm", ""),
                ev.get("tier_mix", ""), ev.get("db_mix", ""),
                "{:,}".format(x0), "{:,}".format(x1)),
             ha="center", va="center", fontsize=6.6, color="#555555")


# ---------------------------------------------------------------------- main
def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--table1", required=True)
    ap.add_argument("--table3", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--spectra", nargs="+", required=True)
    ap.add_argument("--cache", required=True)
    ap.add_argument("--predict", nargs="*", default=[])
    ap.add_argument("--iso-out", required=True)
    ap.add_argument("--igv-out", required=True)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--only-iso", action="store_true")
    ap.add_argument("--only-igv", action="store_true")
    ap.add_argument("--fna", default="",
                    help="genomic.fna（可选，用于在 IGV 页标注注释起始密码子）")
    ap.add_argument("--tx-gtf", default="",
                    help="StringTie merged GTF（转录本结构，用于 RNA-seq 轨道）")
    ap.add_argument("--rna-gtf", nargs="*", default=[],
                    help="StringTie 单样本 GTF（文件/目录），用于转录本 TPM")
    ap.add_argument("--coverage-bedgraph", default="",
                    help="逐碱基深度 bedgraph（bedtools genomecov -bga 产出，可选）")
    args = ap.parse_args(argv)

    t1 = load_table(args.table1)
    t3 = load_table(args.table3)
    with open(args.cache, "rb") as fh:
        cache = pickle.load(fh)

    wanted = {iln(r["peptide"]) for r in t1}
    psms = best_psms(args.spectra, wanted)
    print("peptides %d | with PSM %d" % (len(t1), len(psms)))

    pred = {}
    for p in args.predict:
        try:
            pred.update(read_predictions(p))
        except OSError as exc:
            print("predict skip: %s (%s)" % (p, exc))
    print("prediction keys %d" % len(pred))

    # ---- isoPSV-style pages
    if not args.only_igv:
        rows = sorted(t1, key=lambda r: (r.get("event_type", ""),
                                         r.get("event_sub", ""),
                                         r.get("orf", ""),
                                         int(r.get("pep_nt_start") or 0)))
        with PdfPages(args.iso_out) as pdf:
            fig = plt.figure(figsize=(8.27, 5.85))
            for i, r in enumerate(rows, 1):
                if args.limit and i > args.limit:
                    break
                key = iln(r["peptide"])
                psm = psms.get(key)
                if psm is None:
                    continue
                stem, scan = split_title(psm.get("File_Name", ""))
                peaks = cache.get((stem, scan)) or cache.get((stem, int(scan) if scan.isdigit() else scan)) or []
                r["_type"] = ""
                r["_judgment"] = JUDGMENT.get(str(r.get("tier", "")).strip(),
                                              "only one protease evidence")
                psm = dict(psm)
                psm["_basename"] = psm.get("File_Name", "").split("/")[-1]
                try:
                    chg = int(float(psm.get("Charge") or 0))
                except ValueError:
                    chg = 0
                iso_page(fig, r, psm, peaks, pred.get((r["peptide"], chg)))
                pdf.savefig(fig, dpi=150)
                if i % 50 == 0:
                    print("  iso %d/%d" % (i, len(rows)), flush=True)
            plt.close(fig)
        print("wrote %s" % args.iso_out)

    # ---- IGV-style pages
    if not args.only_iso:
        cds = parse_gtf(args.gtf)
        print("CDS %d" % len(cds))
        by_orf = {}
        for r in t1:
            if not r.get("mapped") or not r.get("pep_nt_start"):
                continue
            by_orf.setdefault(r["orf"], []).append(r)
        events = sorted(t3, key=lambda r: r.get("event_id", ""))

        # 同座位（同 contig/同链/同读框/重叠 ≥80%）的其它事件 ORF：它们的肽在页面上淡色显示，
        # 避免"同一物理基因被拆成两页、各页只看得见一部分肽"的错觉。
        # 必须加同读框判定——不同 frame 的 ORF 是另一个基因，不能当别名。
        def _same_locus(a, b, frac=0.8):
            if a.get("seqid") != b.get("seqid") or a.get("strand") != b.get("strand"):
                return False
            s1, e1 = int(a["orf_start"]), int(a["orf_end"])
            s2, e2 = int(b["orf_start"]), int(b["orf_end"])
            lo, hi = max(s1, s2), min(e1, e2)
            if hi < lo:
                return False
            if a.get("strand") == "-":
                if (e1 - e2) % 3 != 0:
                    return False
            elif (s1 - s2) % 3 != 0:
                return False
            return (hi - lo + 1) / max(1, min(e1 - s1 + 1, e2 - s2 + 1)) >= frac

        alias = {}
        for i, ev in enumerate(events):
            peas = []
            for j, o in enumerate(events):
                if i != j and _same_locus(ev, o):
                    peas.extend(by_orf.get(o["orf"], []))
            alias[ev["orf"]] = peas

        txs, models = [], []
        if args.rna_gtf:
            from orfscape.core.rnaabund import assembled_models, load_transcripts
            models = assembled_models(args.rna_gtf, log=print)
            if args.tx_gtf:                      # 显式给 merged GTF 时才用它
                from orfscape.core.rnaabund import load_samples
                txs = load_transcripts(args.tx_gtf, load_samples(args.rna_gtf))
            else:
                txs = models
            print("转录本行来源：%s（%d 条）"
                  % ("merged.gtf" if args.tx_gtf else "单样本组装", len(txs)))

        depth = {}
        if args.coverage_bedgraph and os.path.exists(args.coverage_bedgraph):
            for line in open(args.coverage_bedgraph, encoding="utf-8", errors="replace"):
                c = line.rstrip(chr(10)).split(chr(9))
                if len(c) < 4:
                    continue
                try:
                    depth.setdefault(c[0], []).append((int(c[1]) + 1, int(c[2]), float(c[3])))
                except ValueError:
                    continue
            for k in depth:
                depth[k].sort()
            print("逐碱基深度：%d 条 contig，%d 段"
                  % (len(depth), sum(len(v) for v in depth.values())))

        genome = {}
        if args.fna and os.path.exists(args.fna):
            name, buf = None, []
            with open(args.fna, encoding="utf-8", errors="replace") as fh:
                for line in fh:
                    if line.startswith(">"):
                        if name:
                            genome[name] = "".join(buf)
                        name, buf = line[1:].split()[0], []
                    else:
                        buf.append(line.strip())
            if name:
                genome[name] = "".join(buf)
            print("genome %d contigs" % len(genome))

        with PdfPages(args.igv_out) as pdf:
            fig = plt.figure(figsize=(11.69, 6.2))
            for i, ev in enumerate(events, 1):
                if args.limit and i > args.limit:
                    break
                os_, oe = int(ev["orf_start"]), int(ev["orf_end"])
                span = oe - os_
                hw = max(int(span * 0.25), 1500)
                win = (max(1, os_ - hw), oe + hw)
                ev["_page"] = i
                try:
                    from igv_v2 import igv_page_v2
                    _sq = genome.get(ev.get("seqid"))
                    igv_page_v2(fig, ev, cds, by_orf.get(ev["orf"], []), win,
                                alt_peps=alias.get(ev["orf"], []),
                                seq=_sq, clen=(len(_sq) if _sq else None),
                                txs=[t for t in txs if t.get("contig") == ev.get("seqid")],
                                models=[m for m in models if m.get("contig") == ev.get("seqid")],
                                depth=depth.get(ev.get("seqid")))
                except ImportError:
                    igv_page(fig, ev, cds, by_orf.get(ev["orf"], []), win)
                pdf.savefig(fig, dpi=150)
                if i % 50 == 0:
                    print("  igv %d/%d" % (i, len(events)), flush=True)
            plt.close(fig)
        print("wrote %s" % args.igv_out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
