# -*- coding: utf-8 -*-
"""输入解析：raw 递归发现、文件名分组推断（馏分制）、镜像组配对、分组表读写。"""
import csv
import re
from pathlib import Path

from ..config import ENZYMES, MIRROR_PAIRS, RAW_PATTERNS

GROUP_FIELDS = ["raw_file", "species_code", "species", "sample_type", "condition",
                "enzyme", "enzyme_name", "cleavage", "fraction",
                "mirror_group_id", "pattern", "group_dir"]


def safe_segment(name):
    """校验单段路径名（物种代码 / 文件名）：禁止分隔符与父目录跳转。"""
    s = str(name or "").strip()
    if (not s or s in (".", "..") or any(c in s for c in ("/", "\\"))
            or Path(s).is_absolute()):
        raise ValueError("unsafe path segment: %r" % (name,))
    return s


def project_species(root):
    """项目里的物种代码：00_input/species.tsv 优先，其次 02_db/05_novel 目录。

    下游（报告、图、收口、P6b）都从这里取物种表，不再写死项目专属代码。
    返回值已过 :func:`safe_segment`，可直接用于拼接路径。
    """
    root = Path(root)
    codes = []

    def add(code):
        try:
            code = safe_segment(code)
        except ValueError:
            return
        if code and code not in codes:
            codes.append(code)

    tsv = root / "00_input" / "species.tsv"
    if tsv.is_file():
        with tsv.open(encoding="utf-8-sig", newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                add((r.get("species_code") or "").strip())
    if not codes:
        for sub in ("02_db", "05_novel", "04_search_annot"):
            d = root / sub
            if not d.is_dir():
                continue
            for p in sorted(d.iterdir()):
                if p.is_dir() and not p.name.startswith((".", "_")):
                    add(p.name)
            if codes:
                break
    return codes


def discover_raws(root):
    """递归发现全部 .raw 文件，返回按路径排序的 Path 列表。"""
    root = Path(root)
    return sorted(p for p in root.rglob("*.raw") if p.is_file())


def infer_one(fname):
    """对单个文件名尝试模式库，返回分组字段 dict 或 None。"""
    for pat in RAW_PATTERNS:
        m = pat["rx"].search(fname)
        if not m:
            continue
        d = m.groupdict()
        sp_code = (pat.get("species")
                   or pat.get("map_species", {}).get(d.get("sp"), d.get("sp")))
        row = {
            "species_code": sp_code,
            "sample_type": pat.get("sample_type") or pat.get("map_sample", {}).get(d.get("st"), d.get("st", "")),
            "enzyme": pat.get("map_enz", {}).get(d.get("enzs"), d.get("enz", d.get("enzs", ""))),
            "replicate": d.get("rep", ""),
            "pattern": pat["name"],
        }
        row["condition"] = row["sample_type"]
        if row["enzyme"] not in ENZYMES:
            row["pattern"] += f"[UNKNOWN_ENZYME:{row['enzyme']}]"
        return row
    return None


def _acq_key(fname):
    """馏分排序键：优先文件名首日期串（XP_ZY_YYYYMMDD_...），回退空。"""
    m = re.search(r"(\d{8})_(\d{8})", fname)
    return m.group(1) if m else ""


def build_groups(raw_paths):
    """批量推断并配对镜像组（馏分制语义，v1.2）。

    语义（用户 2026-10-31 修订）：
    - 每个 (species_code, sample_type, enzyme) 组 = 一个样品组；组内每个 raw = 一个馏分
      （fraction，按文件名采集日期+文件名排序编号 F1..Fn）。技术重跑等同额外馏分。
    - 镜像配对在样品级：同 (species_code, sample_type) 下存在 MIRROR_PAIRS 酶对
      → 分配 mirror_group_id（MG001 起）。
    - 预拆分检测：文件名以 _<CV>.raw 结尾（如 _45.raw/_65.raw，FreeStyle 手工拆分产物）
      标记 pre_split_cv；其对应的多 CV 原 raw（同 stem 去后缀）若同时在场 → 标记
      superseded（后续流程跳过，避免双搜）。

    返回 (rows, unmatched, mirror_mode)。
    """
    rows, unmatched = [], []
    for p in raw_paths:
        g = infer_one(p.name)
        if g is None:
            unmatched.append(p)
            continue
        g.update({
            "raw_file": p.name,
            "path": str(p),
            "group_dir": str(p.parent),
        })
        # 预拆分电压后缀（_45.raw / _-45.raw / _65.raw）。
        # 守卫：|CV| 须在 20–100（FAIMS 合理区间），避免 Mixture_1..6 / M1..M7 这类
        # 馏分编号误判为电压后缀。
        g["pre_split_cv"] = ""
        m_cv = re.search(r"_(-?\d{1,3})\.raw$", p.name)
        if m_cv and 20 <= abs(int(m_cv.group(1))) <= 100:
            g["pre_split_cv"] = m_cv.group(1)
        rows.append(g)

    # 预拆分子 ↔ 原 raw 对应：stem 去掉 _CV 后缀即原文件 stem
    pre_split_stems = {re.sub(r"_(-?\d{1,3})$", "", r["raw_file"][:-4]) for r in rows if r["pre_split_cv"]}
    for r in rows:
        stem = r["raw_file"][:-4]
        r["superseded"] = "yes" if (not r["pre_split_cv"] and stem in pre_split_stems) else ""

    # 馏分编号：组内按 (采集日期, 文件名) 排序 F1..Fn
    groups = {}
    for r in rows:
        groups.setdefault((r["species_code"], r["sample_type"], r["enzyme"]), []).append(r)
    for key in sorted(groups):
        members = sorted(groups[key], key=lambda m: (_acq_key(m["raw_file"]), m["raw_file"]))
        for i, m in enumerate(members, 1):
            m["fraction"] = f"F{i:02d}"
            m["n_fraction_in_group"] = len(members)

    # 样品级镜像配对
    sample_enz = {}
    for r in rows:
        sample_enz.setdefault((r["species_code"], r["sample_type"]), set()).add(r["enzyme"])
    mg_no = 0
    mirror_mode = False
    mg_map = {}
    for key in sorted(sample_enz):
        enz_fs = frozenset(sample_enz[key])
        if any(enz_fs == mp for mp in MIRROR_PAIRS):
            mg_no += 1
            mg_map[key] = f"MG{mg_no:03d}"
            mirror_mode = True
    for r in rows:
        r["mirror_group_id"] = mg_map.get((r["species_code"], r["sample_type"]), "")

    # 物种码 → 酶名/切点信息补全
    for r in rows:
        e = ENZYMES.get(r["enzyme"], {})
        r["enzyme_name"] = e.get("name", r["enzyme"])
        r["cleavage"] = e.get("cleavage", "")
    return rows, unmatched, mirror_mode


def save_groups_csv(rows, out_csv):
    """分组表落盘（用户可手改后再喂回）。"""
    with open(out_csv, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=GROUP_FIELDS + ["path"], extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in GROUP_FIELDS + ["path"]})


def load_groups_csv(in_csv):
    rows = []
    with open(in_csv, newline="", encoding="utf-8-sig") as fh:
        for r in csv.DictReader(fh):
            rows.append(r)
    return rows
