# -*- coding: utf-8 -*-
"""预测器基线：用已注释肽的谱图比较 AlphaPeptDeep / pDeep2 / pDeep3 / ms2pip。

口径
----
* 样本：E1（db_annot）target PSM，每 (肽, 电荷) 取 Final_Score 最优；
* 观测：`03_mgf` 里同 raw/scan 的峰；
* 预测：四个预测器各跑一遍；**统一只取 1+ b/y**（ms2pip v4 CLI 只给 1+，
  四家同口径才可比），20 ppm 内取最强观测峰，sqrt 余弦（`spectra_pdfs.theory_cosine`）；
* 结论：**注释肽上的余弦中位数最高者**即"当下最优预测器"，用于评新肽。

输出（默认写 `04a_qc/`）：
  `predictor_baseline.tsv`  逐 PSM：sid/peptide/charge/mod/cos_<pred>…
  `predictor_summary.tsv`   每个预测器 n / median / IQR / chosen
  `predictor_baseline.png`  箱线 + 散点
"""
import csv
import io
import json
import os
import subprocess
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from .. import localpaths

# 外部工具（可用 ORFSCAPE_* 覆盖；isoPSV 走 localpaths，兜底用通用目录）
ISO = Path(localpaths.get("isopsv", r"D:\isoPSV_multiengine",
                          os.path.join(os.path.expanduser("~"), "isoPSV_multiengine")))
PFIND_PY = os.environ.get("ORFSCAPE_ALPHA_PY", "D:/Bioinformatics/Miniconda/python.exe")
PY37 = os.environ.get("ORFSCAPE_PDEEP_PY", "D:/Bioinformatics/pDeep2/py37/python.exe")
PDEEP2_DIR = Path(os.environ.get("ORFSCAPE_PDEEP2_DIR",
                                 "D:/Bioinformatics/pDeep2/pDeep/pDeep2"))
MS2PIP = Path(os.environ.get("ORFSCAPE_MS2PIP",
                             "D:/Bioinformatics/ms2pip/venv/Scripts/ms2pip.exe"))

BLUE, ORANGE, GREY, RED = "#2F6FB0", "#E07B39", "#8C8C8C", "#B2182B"
plt.rcParams.update({"font.family": "Arial", "font.size": 8,
                     "axes.linewidth": 0.7, "axes.spines.top": False,
                     "axes.spines.right": False, "savefig.dpi": 300})


def _rows(path):
    with open(path, encoding="utf-8-sig", newline="", errors="replace") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_mgf_block(path, wanted):
    """读指定 scan 的峰；返回 {scan: (mz[], ints[])}。"""
    out = {}
    title, buf, in_block = None, [], False
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            s = line.strip()
            if s == "BEGIN IONS":
                in_block, title, buf = True, None, []
                continue
            if not in_block:
                continue
            if s == "END IONS":
                in_block = False
                scan = None
                if title:
                    # TITLE 两种写法都要认：`stem.scan.scan.chg.0.dta`（msconvert）
                    # 与 `stem.scan.scan.chg`（pParse 等）。只认前者会让观测谱恒为 0 条。
                    t = title[:-4] if title.endswith(".dta") else title
                    parts = t.split(".")
                    if len(parts) >= 4 and parts[-4].isdigit():
                        scan = int(parts[-4])
                    elif len(parts) >= 3 and parts[-3].isdigit():
                        scan = int(parts[-3])
                if scan in wanted and buf:
                    arr = np.array(buf, dtype=float)
                    out[scan] = (arr[:, 0], arr[:, 1])
                continue
            if s.startswith("TITLE="):
                title = s[6:]
            elif s and not s[0].isalpha() and not s.startswith("CHARGE"):
                parts = s.split()
                if len(parts) >= 2:
                    try:
                        buf.append((float(parts[0]), float(parts[1])))
                    except ValueError:
                        pass
    return out


def pick_sample(workdir, n, log=print):
    """E1 注释 PSM：每 (肽, 电荷) 取 Final_Score 最优。"""
    best = {}
    for spec in Path(workdir, "04_search_annot").glob("*/" + "*/pFind-Filtered.spectra"):
        for r in _rows(spec):
            if str(r.get("Target/Decoy", "")).lower().startswith("dec"):
                continue
            pep = (r.get("Sequence") or "").upper().replace("I", "L")
            if not pep:
                continue
            try:
                score = float(r.get("Final_Score") or 0)
            except ValueError:
                continue
            key = (pep, r.get("Charge"))
            cur = best.get(key)
            if cur is None or score < cur["score"]:
                fname = str(r.get("File_Name") or "")
                stem = fname.rstrip(".dta").split(".")[0]
                best[key] = dict(sid="%s.%s" % (stem, r.get("Scan_No")), pep=pep,
                                 mod=(r.get("Modification") or "").strip(),
                                 charge=int(float(r.get("Charge") or 2)),
                                 stem=stem, scan=int(r.get("Scan_No")), score=score)
    items = sorted(best.values(), key=lambda d: d["sid"])
    if n and len(items) > n:
        step = len(items) / float(n)
        items = [items[int(i * step)] for i in range(n)]
    log("  注释肽样本：%d（候选 %d）" % (len(items), len(best)))
    return items


def convert_mod(mod):
    """pFind '2,Oxidation[M];' → peprec '2|Oxidation|13|Carbamidomethyl'。"""
    out = []
    for tok in (mod or "").split(";"):
        tok = tok.strip()
        if not tok or "," not in tok:
            continue
        pos, name = tok.split(",", 1)
        base = name.split("[")[0].strip()
        if base:
            out.append("%s|%s" % (pos, base))
    return "|".join(out)


def _run(cmd, cwd=None, env=None, timeout=3600):
    return subprocess.run(cmd, cwd=cwd, env=env, capture_output=True, text=True,
                          encoding="utf-8", errors="replace", timeout=timeout)


def run_alpha(items, tmp):
    inp, out = tmp / "alpha_in.tsv", tmp / "alpha_out.tsv"
    lines = ["spec_id\tpeptide\tmodinfo\tcharge"]
    lines += ["%s\t%s\t%s\t%d" % (d["sid"], d["pep"], d["mod"], d["charge"])
              for d in items]
    inp.write_text("\n".join(lines) + "\n", encoding="utf-8")
    r = _run([PFIND_PY, str(ISO / "code/predictors/alphapeptdeep_driver.py"),
              "--in", str(inp), "--out", str(out), "--instrument", "Lumos",
              "--ce", "30"])
    if not out.is_file():
        return None, (r.stderr or r.stdout or "")[-300:]
    return {x["spec_id"]: x for x in _rows(out)}, ""


def run_pdeep3(items, tmp):
    inp, out = tmp / "pdeep3_in.tsv", tmp / "pdeep3_out.tsv"
    lines = ["spec_id\tpeptide\tmodinfo\tcharge"]
    lines += ["%s\t%s\t%s\t%d" % (d["sid"], d["pep"], d["mod"], d["charge"])
              for d in items]
    inp.write_text("\n".join(lines) + "\n", encoding="utf-8")
    env = dict(os.environ)
    env.setdefault("PDEEP3_ROOT", "D:/Bioinformatics/pDeep3")
    r = _run([PY37, str(ISO / "code/predictors/pdeep3_driver.py"),
              "--in", str(inp), "--out", str(out), "--instrument", "Lumos",
              "--ce", "30"], env=env)
    if not out.is_file():
        return None, (r.stderr or r.stdout or "")[-300:]
    return {x["spec_id"]: x for x in _rows(out)}, ""


def run_pdeep2(items, tmp):
    inp, out = tmp / "pdeep2_in.txt", tmp / "pdeep2_out.txt"
    lines = ["%s\t%s\t%d" % (d["pep"], d["mod"], d["charge"]) for d in items]
    inp.write_text("\n".join(lines) + "\n", encoding="utf-8")
    r = _run([PY37, "predict.py", "--nce", "0.30", "--instrument", "Lumos",
              "--input", str(inp), "--output", str(out)], cwd=str(PDEEP2_DIR))
    if not out.is_file():
        return None, (r.stderr or r.stdout or "")[-300:]
    res, cur = {}, None
    with open(out, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("TITLE="):
                parts = line[6:].strip().split("|")
                if len(parts) >= 2:
                    cur = (parts[0], parts[-1])
                    res.setdefault(cur, {})
                else:
                    cur = None
            elif cur and line.startswith(("b+1=", "b+2=", "y+1=", "y+2=")):
                k, v = line.split("=", 1)
                res[cur][k] = v.strip()
            elif line.startswith("END IONS"):
                cur = None
    return res, ""


def run_ms2pip(items, tmp):
    inp, out = tmp / "ms2pip_in.tsv", tmp / "ms2pip_out.tsv"
    lines = ["spec_id\tpeptide\tmodifications\tcharge"]
    lines += ["%s\t%s\t%s\t%d" % (d["sid"], d["pep"], convert_mod(d["mod"]),
                                 d["charge"]) for d in items]
    inp.write_text("\n".join(lines) + "\n", encoding="utf-8")
    r = _run([str(MS2PIP), "predict-batch", str(inp), "-t", "peprec", "-f", "tsv",
              "-o", str(tmp / "ms2pip_out"), "--model", "HCD"])
    if not out.is_file():
        return None, (r.stderr or r.stdout or "")[-300:]
    series = {}
    for x in _rows(out):
        try:
            i = int(x["psm_index"])
        except (TypeError, ValueError):
            continue
        if i >= len(items):
            continue
        d = series.setdefault(items[i]["sid"], {"b+1": {}, "y+1": {}})
        key = "b+1" if x["ion_type"] == "b" else "y+1"
        try:
            d[key][int(x["ion_number"])] = 2.0 ** float(x["predicted"])
        except (ValueError, OverflowError):
            continue
    out2 = {}
    for sid, d in series.items():
        obj = {}
        for key, nums in d.items():
            n = max(nums) if nums else 0
            if n:
                arr = [0.0] * n
                for k, v in nums.items():
                    arr[k - 1] = v
                obj[key] = arr
        out2[sid] = obj
    return out2, ""


def cos_for(seq, mod, peaks, series, theory):
    if not series or not peaks:
        return None
    mz, ints = peaks
    sticks = theory.pred_sticks(seq, mod, {"b+1": series.get("b+1"),
                                           "y+1": series.get("y+1")})
    if not sticks:
        return None
    return theory.theory_cosine(mz, ints, sticks, 20.0)


def _series_from_tsv(raw):
    ser = {}
    if not raw:
        return None
    if raw.get("b1"):
        ser["b+1"] = [float(x) for x in raw["b1"].split(",")]
    if raw.get("y1"):
        ser["y+1"] = [float(x) for x in raw["y1"].split(",")]
    return ser or None


def _series_from_pdeep2(raw):
    ser = {}
    if not raw:
        return None
    for k in ("b+1", "y+1"):
        if raw.get(k):
            ser[k] = [float(x) for x in raw[k].split(",")]
    return ser or None


def run(workdir, out_dir=None, n=600, log=print):
    """跑一遍基线；返回 summary dict（含 chosen）。"""
    import importlib
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
    theory = importlib.import_module("spectra_pdfs")

    wd = Path(workdir).resolve()
    out = Path(out_dir).resolve() if out_dir else wd / "04a_qc"
    out.mkdir(parents=True, exist_ok=True)
    tmp = out / "_inputs"
    tmp.mkdir(exist_ok=True)

    items = pick_sample(wd, n, log=log)
    by_stem = {}
    for d in items:
        by_stem.setdefault(d["stem"], set()).add(d["scan"])
    peaks = {}
    for stem, scans in by_stem.items():
        mgf = Path(wd, "03_mgf", stem + ".mgf")
        if mgf.is_file():
            peaks[stem] = read_mgf_block(str(mgf), scans)
    log("  观测谱：%d/%d" % (sum(len(v) for v in peaks.values()), len(items)))

    runners = [("alphapeptdeep", run_alpha), ("pdeep3", run_pdeep3),
               ("pdeep2", run_pdeep2), ("ms2pip", run_ms2pip)]
    results, errors = {}, {}
    for name, fn in runners:
        try:
            res, err = fn(items, tmp)
        except Exception as exc:                       # noqa: BLE001
            res, err = None, "%s: %s" % (type(exc).__name__, exc)
        if res is None:
            errors[name] = err
            log("  %-13s 跳过：%s" % (name, err[:100].replace("\n", " ")))
            continue
        results[name] = res
        log("  %-13s 预测 %d 条" % (name, len(res)))

    table = []
    for d in items:
        spec = peaks.get(d["stem"], {}).get(d["scan"])
        row = dict(sid=d["sid"], peptide=d["pep"], charge=d["charge"], mod=d["mod"])
        for name in results:
            if name in ("alphapeptdeep", "pdeep3"):
                ser = _series_from_tsv(results[name].get(d["sid"]))
            elif name == "pdeep2":
                ser = _series_from_pdeep2(results[name].get((d["pep"], str(d["charge"]))))
            else:
                ser = results[name].get(d["sid"])
            row["cos_" + name] = cos_for(d["pep"], d["mod"], spec, ser, theory)
        table.append(row)

    fields = ["sid", "peptide", "charge", "mod"] + ["cos_" + n for n in results]
    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=fields, delimiter="\t")
    w.writeheader()
    for r in table:
        w.writerow({k: ("" if r.get(k) is None else
                        ("%.4f" % r[k] if isinstance(r.get(k), float) else r.get(k)))
                    for k in fields})
    (out / "predictor_baseline.tsv").write_text(buf.getvalue(), encoding="utf-8")

    summary = {}
    for name in results:
        vals = [r["cos_" + name] for r in table if r.get("cos_" + name) is not None]
        if not vals:
            continue
        arr = np.array(vals)
        summary[name] = dict(n=int(arr.size), median=float(np.median(arr)),
                             q1=float(np.percentile(arr, 25)),
                             q3=float(np.percentile(arr, 75)))
    chosen = max(summary, key=lambda k: summary[k]["median"]) if summary else ""
    lines = ["predictor\tn\tmedian\tq1\tq3\tchosen\terror"]
    for name in ("alphapeptdeep", "pdeep2", "pdeep3", "ms2pip"):
        s = summary.get(name)
        if s:
            lines.append("%s\t%d\t%.4f\t%.4f\t%.4f\t%s\t"
                         % (name, s["n"], s["median"], s["q1"], s["q3"],
                            "1" if name == chosen else "0"))
        else:
            lines.append("%s\t0\t\t\t\t0\t%s" % (name, errors.get(name, "")[:120]))
    (out / "predictor_summary.tsv").write_text("\n".join(lines) + "\n",
                                               encoding="utf-8")
    log("  注释肽基线：%s（chosen=%s）" % (json.dumps(summary, ensure_ascii=False), chosen))

    if summary:
        fig, ax = plt.subplots(figsize=(5.4, 3.0))
        order = [n for n in ("alphapeptdeep", "pdeep2", "pdeep3", "ms2pip")
                 if n in summary]
        data = [[r["cos_" + n] for r in table if r.get("cos_" + n) is not None]
                for n in order]
        bp = ax.boxplot(data, widths=0.55, patch_artist=True,
                        medianprops=dict(color="k", lw=1.0),
                        flierprops=dict(marker="", ls="none"))
        for patch, c in zip(bp["boxes"], (BLUE, ORANGE, GREY, RED)):
            patch.set_facecolor(c)
            patch.set_edgecolor("#333333")
            patch.set_linewidth(0.6)
        rng = np.random.default_rng(3)
        for i, vals in enumerate(data):
            xs = i + 1 + rng.uniform(-0.16, 0.16, len(vals))
            ax.scatter(xs, vals, s=1.0, color="#555555", alpha=0.30, lw=0)
        ax.set_xticks(range(1, len(order) + 1))
        ax.set_xticklabels(["%s\nn=%d" % (n, summary[n]["n"]) for n in order],
                           fontsize=7)
        ax.set_ylabel("cosine (1+ b/y)")
        ax.set_title("Annotated-peptide baseline", loc="left", fontsize=9,
                     fontweight="bold")
        fig.tight_layout()
        png = out / "predictor_baseline.png"
        fig.savefig(str(png), bbox_inches="tight", facecolor="white")
        try:
            fig.savefig(str(out / "predictor_baseline.pdf"), bbox_inches="tight",
                        facecolor="white")
        except OSError:
            pass
        plt.close(fig)
    return dict(predictors=summary, chosen=chosen, errors=errors, n=len(items))
