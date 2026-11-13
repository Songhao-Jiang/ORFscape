# -*- coding: utf-8 -*-
"""碎片离子质量误差（ppm）自算。

背景：pFind 谱表 `Avg.Frag.Mass.Shift` 列的实测均值与 summary
`Fragment_Mean` 完全相同（-1.26），但 -1.26 Da 作为碎片误差物理上不成立
（真实碎片误差为 ppm 级）；独立复算显示该列量级为 mDa/ppm 而非 Da。为使
E1Q4 的母离子与子离子两面板单位统一为 ppm，本模块从 03_mgf 原始峰与 b/y
理论离子重算。

安全约定：MGF 只经 03_mgf 目录的 scandir 白名单访问，文件描述符由 os.open
取得，不做数据派生的路径拼接；本模块不写任何文件。
"""
import bisect
import os

RES_MASS = dict(A=71.03711, R=156.10111, N=114.04293, D=115.02694,
                C=103.00919, E=129.04259, Q=128.05858, G=57.02146,
                H=137.05891, I=113.08406, L=113.08406, K=128.09496,
                M=131.04049, F=147.06841, P=97.05276, S=87.03203,
                T=101.04768, W=186.07931, Y=163.06333, V=99.06841)
PROTON_M, H2O_M = 1.007276, 18.010565
MOD_DELTA = (("Carbamidomethyl", 57.021464), ("Oxidation", 15.994915),
             ("Acetyl", 42.010565))
TOL_PPM = 20.0
MIN_IONS = 4


def parse_mods(mod):
    """pFind Modification → (N端增量, {位点1based: 增量})。"""
    nterm, pos = 0.0, {}
    for part in filter(None, (mod or "").split(";")):
        f = part.split(",")
        if len(f) < 2:
            continue
        try:
            p = int(f[0])
        except ValueError:
            continue
        for name, dv in MOD_DELTA:
            if name in f[1]:
                if name == "Acetyl" and p <= 1:
                    nterm += dv
                else:
                    pos[p] = pos.get(p, 0) + dv
                break
    return nterm, pos


def theo_ions(seq, nterm, mods):
    """b/y 单双电荷理论 m/z 列表。"""
    rm = [RES_MASS.get(a, 0.0) + mods.get(i + 1, 0.0)
          for i, a in enumerate(seq)]
    if rm:
        rm[0] += nterm
    out = []
    pre = 0.0
    for i in range(len(seq) - 1):
        pre += rm[i]
        for z in (1, 2):
            out.append((pre + z * PROTON_M) / z)
    suf = 0.0
    for i in range(len(seq) - 1, 0, -1):
        suf += rm[i]
        base = suf + H2O_M
        for z in (1, 2):
            out.append((base + z * PROTON_M) / z)
    return out


def _ppm_of(peaks, seq, mod):
    """单谱：理论离子与实测峰匹配后的 ppm 均值（匹配不足返回 None）。"""
    mzs = sorted(p[0] for p in peaks)
    if not mzs:
        return None
    nm, pm = parse_mods(mod)
    errs = []
    for t in theo_ions(seq, nm, pm):
        tol = t * TOL_PPM * 1e-6
        j = bisect.bisect_left(mzs, t - tol)
        best = None
        while j < len(mzs) and mzs[j] <= t + tol:
            d = abs(mzs[j] - t)
            if best is None or d < best:
                best = mzs[j]
            j += 1
        if best is not None:
            errs.append((best - t) / t * 1e6)
    if len(errs) < MIN_IONS:
        return None
    return sum(errs) / len(errs), len(errs)


def _norm_title(title):
    """`stem.100.100.2` / `stem.100.100.2.0.dta` → 统一四段形态（pParse 差异）。"""
    t = (title or "").strip()
    if t.lower().endswith(".dta"):
        t = t[:-4]
    parts = t.split(".")
    if len(parts) >= 5 and parts[-1] == "0":
        t = ".".join(parts[:-1])
    return t


def _scan_handle(fh, wanted, sink):
    """流式解析一个已打开的 MGF，命中 wanted 的谱写入 sink。"""
    title, peaks = None, None
    for line in fh:
        s = line.rstrip("\n")
        if s == "BEGIN IONS":
            title, peaks = None, []
        elif s.startswith("TITLE=") and peaks is not None:
            title = s[6:].strip()
        elif s == "END IONS":
            if title and peaks:
                parts = _norm_title(title).rsplit(".", 3)
                hit = (wanted.get(".".join(parts[1:]))
                       if len(parts) == 4 else None)
                if hit:
                    res = _ppm_of(peaks, hit[1], hit[2])
                    if res:
                        sink[hit[0]] = res
            title, peaks = None, None
        elif peaks is not None and s and s[0].isdigit():
            f = s.split()
            if len(f) >= 2:
                try:
                    peaks.append((float(f[0]), float(f[1])))
                except ValueError:
                    pass


def compute(root, psm_meta, log=print):
    """psm_meta: {File_Name: (seq, modification)} → {File_Name: (ppm, n_ions)}。"""
    by_stem = {}
    for fn, meta in psm_meta.items():
        parts = _norm_title(fn).rsplit(".", 3)
        if len(parts) != 4:
            continue
        by_stem.setdefault(parts[0], {})[".".join(parts[1:])] = (fn,) + tuple(meta)
    sink = {}
    scan_root = os.path.join(root, "03_mgf")
    if os.path.isdir(scan_root):
        with os.scandir(scan_root) as it:
            for entry in it:
                if not (entry.is_file()
                        and entry.name.lower().endswith(".mgf")):
                    continue
                wanted = by_stem.get(entry.name[:-4])
                if not wanted:
                    continue
                fd = os.open(entry.path, os.O_RDONLY)
                with os.fdopen(fd, "r", encoding="utf-8",
                               errors="replace") as fh:
                    _scan_handle(fh, wanted, sink)
    log(f"  碎片 ppm 自算完成：{len(sink):,} / {len(psm_meta):,} 谱")
    return sink
