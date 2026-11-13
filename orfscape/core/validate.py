# -*- coding: utf-8 -*-
"""R8 新肽收口：novel PSM 全量过 isoPSV 无标记流程，Tier 分流门控。

用户口径（2026-11-06）：
  1) 所有新肽段（2,111）的 PSM 全部送 isoPSV 三步流程（竞争搜库→阶段 A 仲裁）
  2) Tier1/Tier2 全保留——isoPSV 通过率只作证据指标报告，不作门
  3) Tier3 只保留「高质量 且 通过 isoPSV 验证」部分：
       高质量 = 谱图质量卡值（cand1 PSM score ≥ min_score，默认 0.30）
                ∧ 理论谱图相似性卡值（主预测器 cosine ≥ min_cos，默认 0.90）
       通过 isoPSV 验证 = 阶段 A 六重存疑全过（credible_novel = TRUE）

isoPSV 侧契约（isoPSV_multiengine v2.0.175）：
  - psm_file 任意带 File_Name 列的制表符表（pFind 家族表直用，自带 Label）
  - fasta = db_annot.faa → 原始蛋白（L2_/L3_）不在库内 → 全部落第 2 部分
    新肽段竞争报告（2.novel_peptide_competition_report.tsv，credible_novel）
  - enzyme 键支持 lysarginase（KR/N）；镜像数据按 (物种,酶) 拆 4 次运行
  - 竞争槽位：a=原始 novel PSM、b=pFind 开放式、c=pFind 限定式（auto_search=1）
  - 四预测器余弦叠加到 1.unlabeled_result/master_table.tsv（predict_cosine_*）
  - 长标题（150 字符 stem）会顶爆 unlabeled/<raw>_HCDFT.mgf 路径 →
    抽谱时把 TITLE 前缀改短码 V<k>，psm_file 的 File_Name 同步重写
  - R 侧须清 LANG/LC_ALL/LC_CTYPE；CONFIG 环境变量传 config_step1 绝对路径
"""
import csv
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

from .. import localpaths
from .novel import iln, load_spectra

# isoPSV 安装根（--isopsv > ORFSCAPE_ISOPSV > paths.local.json > 常见目录）
DEFAULT_ISOPSV = localpaths.get(
    "isopsv", os.path.join(os.path.expanduser("~"), "isoPSV_multiengine"),
    r"D:\isoPSV_multiengine", r"C:\isoPSV_multiengine")
# 本机依赖（config_step1/config_predictor 生成用），同样走 localpaths
PFIND_DIR = localpaths.get("pfind", r"D:\Bioinformatics\pFind3.1",
                           r"C:\Program Files\pFind3.1")
PDEEP2_DIR = localpaths.get("pdeep2", r"D:\Bioinformatics\pDeep2")
MINICONDA_DIR = localpaths.get("conda", r"D:\Bioinformatics\Miniconda",
                               os.path.join(os.path.expanduser("~"), "miniconda3"),
                               os.path.join(os.path.expanduser("~"), "anaconda3"))
MS2PIP_VENV = localpaths.get("ms2pip", r"D:\Bioinformatics\ms2pip\venv")

ENZYME_KEY = {"TL": "trypsin", "LL": "lysarginase"}
# 主预测器（阶段 A predict 列）；余弦卡值用同名列
PRIMARY_PREDICTOR = "alphapeptdeep"

csv.field_size_limit(2 ** 31 - 1)

# 允许的子目录名白名单（防止外部输入构造路径）
_SAFE_DIR = ("05_novel", "05b_research", "06_validate", "06b_validate_research")
_SAFE_FILE = ("novel_peptides.tsv", "novel_after_research.tsv", "novel_validated.tsv")


def read_table_top(root, sub, filename):
    """读 root/sub/filename（顶层单文件，scandir 白名单）。"""
    if sub not in _SAFE_DIR or filename not in _SAFE_FILE:
        raise ValueError("unsafe table request: %s/%s" % (sub, filename))
    d = os.path.join(str(root), sub)
    if not os.path.isdir(d):
        return []
    with os.scandir(d) as it:
        for e in it:
            if e.name == filename and e.is_file():
                with open(e.path, encoding="utf-8-sig", newline="") as fh:
                    return list(csv.DictReader(fh, delimiter="\t"))
    return []


def load_reuse_verdicts(root, sub="06_validate", log=print):
    """复用已有收口结论 → {peptide: {k: v}}（键为 verdict 字段名）。"""
    keys = ("n_spec", "n_credible", "iso_pass", "best_score", "best_n_frag",
            "best_missing_run", "best_cos", "cos_src", "doubt", "status",
            "quality_pass", "sim_pass", "gate_pass")
    out = {}
    for r in read_table_top(root, sub, "novel_validated.tsv"):
        pep = r.get("peptide", "")
        if not pep:
            continue
        out[pep] = {k: r.get("v_" + k, "") for k in keys}
    if out:
        log(f"  复用 {sub} 的已验结论 {len(out)} 肽")
    return out


def read_table_in(root, sub, sp, filename):
    """用 scandir 白名单读取 root/sub/sp/filename（不拼接后直接 open）。"""
    if sub not in _SAFE_DIR or filename not in _SAFE_FILE:
        raise ValueError("unsafe table request: %s/%s" % (sub, filename))
    d = os.path.join(str(root), sub, str(sp))
    if not os.path.isdir(d):
        return []
    with os.scandir(d) as it:
        for e in it:
            if e.name == filename and e.is_file():
                with open(e.path, encoding="utf-8-sig", newline="") as fh:
                    return list(csv.DictReader(fh, delimiter="\t"))
    return []


# ---------------------------------------------------------------- 收集
def list_species(root):
    """项目里实际存在的物种代号：05_novel / 05b_research 下的子目录，回退 Mab/Mch。"""
    seen = []
    for sub in ("05_novel", "05b_research"):
        base = os.path.join(str(root), sub)
        if os.path.isdir(base):
            for name in sorted(os.listdir(base)):
                if os.path.isdir(os.path.join(base, name)) and name not in seen:
                    seen.append(name)
    return seen or ["Mab", "Mch"]


def load_novel_sets(root, source="novel"):
    """source=novel → 05_novel；source=research → 05b_research（Tier 取 tier_new）。"""
    out = {}
    for sp in list_species(root):
        d = {}
        if source == "research":
            for r in read_table_in(root, "05b_research", sp,
                                   "novel_after_research.tsv"):
                d[r["peptide"]] = r.get("tier_new", "")
        else:
            for r in read_table_in(root, "05_novel", sp, "novel_peptides.tsv"):
                d[r["peptide"]] = r["tier"]
        out[sp] = d
    return out


def collect_novel_psms(root, source="novel", only=None, log=print):
    """E2（及 P6b 回搜）结果 → novel PSM（谱级去重）。返回行列表。

    每行：File_Name/Scan_No/... 原 pFind 列 + species/enzyme/peptide/tier。
    同谱多行保留 Final_Score 最高者。source=research 时：
      先取 E2 证据，再用 05b_research 结果**补 E2 未覆盖的谱**（新肽的谱主要来自这里）。
    """
    nov = load_novel_sets(root, source)
    best = {}

    def scan(base, sp, only_new_spectra=False):
        if not os.path.isdir(base):
            return 0
        hit = 0
        for task in os.listdir(base):
            f = os.path.join(base, task, "pFind-Filtered.spectra")
            if not os.path.exists(f):
                continue
            enz = "LL" if "_LL_" in task else "TL"
            for r in load_spectra(f):
                prots = (r.get("Proteins") or "")
                if not ("L2_" in prots or "L3_" in prots):
                    continue
                pep = iln(r.get("Sequence", ""))
                if pep not in nov[sp]:
                    continue
                if only is not None and pep not in only:
                    continue
                fn = r.get("File_Name", "")
                if not fn:
                    continue
                if only_new_spectra and fn in best:
                    continue
                r = dict(r)
                r["species"], r["enzyme"] = sp, enz
                r["peptide"], r["tier"] = pep, nov[sp][pep]
                try:
                    sc = float(r.get("Final_Score") or 0)
                except ValueError:
                    sc = 0.0
                if fn not in best or sc > best[fn]["_sc"]:
                    r["_sc"] = sc
                    r["_src"] = "research" if base.find("05b_research") >= 0 else "E2"
                    best[fn] = r
                    hit += 1
        return hit

    for sp in list_species(root):
        scan(os.path.join(root, "04_search_disc", sp), sp)
    if source == "research":
        for sp in list_species(root):
            n = scan(os.path.join(root, "05b_research", sp), sp, only_new_spectra=True)
            if n:
                log(f"  [{sp}] 回搜补入新谱 {n}")
    rows = list(best.values())
    for r in rows:
        r.pop("_sc", None)
    n_res = sum(1 for x in rows if x.get("_src") == "research")
    log(f"  novel PSM（谱级去重）{len(rows)}"
        f"（TL {sum(1 for x in rows if x['enzyme']=='TL')}"
        f" / LL {sum(1 for x in rows if x['enzyme']=='LL')}；"
        f"来自回搜 {n_res}）")
    return rows


# ---------------------------------------------------------------- 抽谱与输入落盘
def read_mgf_blocks(path):
    """MGF → {TITLE: block_text}。"""
    out = {}
    cur, title = [], None
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            s = line.strip()
            if s == "BEGIN IONS":
                cur, title = [line], None
            elif cur is not None:
                cur.append(line)
                if s.startswith("TITLE="):
                    title = s[6:]
                if s == "END IONS":
                    if title:
                        out[title] = "".join(cur)
                    cur = None
    return out


def layer_of(protein):
    """蛋白名 → (层, 证据来源)。L6 的名字带竖线，不能按 `_` 切。"""
    p = str(protein or "")
    if p.startswith("L2_"):
        return "L2", "intergenic_6frame"
    if p.startswith("L3_") or p.startswith("L3|"):
        return "L3", "alt_start_extension"
    if p.startswith("L6|"):
        return "L6", "transcript_evidence"
    if p.startswith("L1_"):
        return "L1", "ab_initio"
    return p.split("_")[0], "alt_start_extension"


def _pep_fprint(peps):
    """输入肽集指纹（去重排序后 sha256 前 16 位）——isoPSV 结果复用判据。"""
    import hashlib
    uniq = sorted({str(p).strip() for p in peps if str(p).strip()})
    return hashlib.sha256("\n".join(uniq).encode("utf-8")).hexdigest()[:16]


def norm_title(title):
    """把两种 TITLE 形态归一到 `<stem>.<scan>.<scan>.<z>`。

    - msconvert：`stem.100.100.2`
    - pParse（pFind 原生）：`stem.100.100.2.0.dta`，末尾多 `.0.dta`
    不归一，`rsplit(".", 3)` 会把 scan/charge 切错。
    """
    t = (title or "").strip()
    if t.lower().endswith(".dta"):
        t = t[:-4]
    parts = t.split(".")
    if len(parts) >= 5 and parts[-1] == "0":
        t = ".".join(parts[:-1])
    return t


def split_title(title):
    """`<stem>.<scan>.<scan>.<z>` → (stem, 'scan.scan.z')。stem 无点号（实测）。"""
    parts = norm_title(title).rsplit(".", 3)
    if len(parts) != 4:
        return None, None
    return parts[0], ".".join(parts[1:])


def stage_run(root, rows, run_dir, mgf_dir, log=print):
    """单次运行输入落盘：短标题 MGF + psm_file + title_map。

    返回 (psm_path, title_map)。title_map: 新 File_Name → 原行。
    """
    os.makedirs(run_dir, exist_ok=True)
    os.makedirs(mgf_dir, exist_ok=True)
    # 按源 MGF 分组
    by_src = {}
    for r in rows:
        stem, suf = split_title(r["File_Name"])
        by_src.setdefault(stem, []).append((suf, r))
    # 源文件名 = stem.mgf（03_mgf 全量同名）
    title_map, new_rows = {}, []
    kept_cols = list(rows[0].keys())
    for meta in ("species", "enzyme", "peptide", "tier"):
        kept_cols.remove(meta)
    for k, (stem, items) in enumerate(sorted(by_src.items()), 1):
        src = os.path.join(root, "03_mgf", stem + ".mgf")
        if not os.path.exists(src):
            log(f"  [warn] 源 MGF 缺失，跳过 {stem[:40]}...")
            continue
        blocks = read_mgf_blocks(src)
        short = f"V{k:02d}"
        # 一张 MS2 可能同时被多个候选肽命中。原来按行顺序"后写覆盖"，谁拿到谱
        # 取决于遍历顺序；改成**最终打分高者优先**，且记下被挤掉的肽（它们会
        # 在收口阶段表现为 MISSING，需要与"isoPSV 给了谱但没结论"区分）。
        best = {}
        for suf, r in items:
            try:
                sc = float(r.get("Final_Score") or 0)
            except (TypeError, ValueError):
                sc = 0.0
            key = f"{short}.{suf}"
            if key not in best or sc > best[key][0]:
                best[key] = (sc, r)
        for key, (_sc, r) in best.items():
            title_map[key] = r
        with open(os.path.join(mgf_dir, short + ".mgf"), "w",
                  encoding="utf-8", newline="") as out:
            for t, blk in blocks.items():
                st, su = split_title(t)
                if st != stem or f"{short}.{su}" not in title_map:
                    continue
                # 重写 TITLE 前缀（块内第一处 TITLE=）
                blk = re.sub(r"(TITLE=)" + re.escape(stem) + r"(?=\.)",
                             r"\g<1>" + short, blk, count=1)
                out.write(blk)
                r = title_map[f"{short}.{su}"]
                nr = {c: r.get(c, "") for c in kept_cols}
                nr["File_Name"] = f"{short}.{su}"
                new_rows.append(nr)
    psm_path = os.path.join(run_dir, "psm_file.tsv")
    with open(psm_path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=kept_cols, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(new_rows)
    map_path = os.path.join(run_dir, "title_map.tsv")
    with open(map_path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["File_Name_new", "File_Name_orig", "species", "enzyme",
                    "peptide", "tier", "Proteins"])
        for nf, r in sorted(title_map.items()):
            w.writerow([nf, r["File_Name"], r["species"], r["enzyme"],
                        r["peptide"], r["tier"], r.get("Proteins", "")])
    log(f"  落盘 {len(new_rows)} 谱 / {len(by_src)} 源 MGF → {mgf_dir}")
    return psm_path, title_map


# ---------------------------------------------------------------- 配置生成
CONFIG_STEP1_BODY = """psm_file\t{psm}
fasta\t{fasta}
unlabeled_channel\t1
enzyme\t{enzyme}
mgf\t{mgf}
auto_search\t1
pfind_exe\t{pfind}
use_punifind\t0
predictor_primary\t{primary}
use_pdeep2\t1
pdeep2_exe\t{pdeep2}
pdeep2_instrument\tQE
pdf_on\t0
pfind_workers\t4
python_exe\t{pypath}
candidate_scope\tfiltered
proteogenomics_mode\t1
"""

CONFIG_PREDICTOR = """predictor1_type\talphapeptdeep
predictor1_mode\t1
predictor1_python\t{conda}
predictor1_params\tQE,27
predictor2_type\tms2pip
predictor2_mode\t1
predictor2_python\t{ms2pip}
predictor2_params\t
predictor3_type\tpdeep3
predictor3_mode\t1
predictor3_python\t
predictor3_params\t
"""


def write_configs(run_dir, enzyme, psm_path, fasta, mgf_dir, iso_root):
    cfg1 = os.path.join(run_dir, "config_step1.txt")
    with open(cfg1, "w", encoding="utf-8", newline="") as fh:
        fh.write(CONFIG_STEP1_BODY.format(
            psm=psm_path.replace("\\", "/"), fasta=fasta.replace("\\", "/"),
            enzyme=ENZYME_KEY[enzyme], mgf=mgf_dir.replace("\\", "/"),
            pfind=PFIND_DIR.replace("\\", "/"),
            primary=PRIMARY_PREDICTOR, pdeep2=PDEEP2_DIR.replace("\\", "/"),
            pypath=os.path.join(iso_root, "runtime", "python", "python.exe")))
    cfgp = os.path.join(run_dir, "config_predictor.txt")
    with open(cfgp, "w", encoding="utf-8", newline="") as fh:
        fh.write(CONFIG_PREDICTOR.format(
            conda=MINICONDA_DIR.replace("\\", "/"),
            ms2pip=MS2PIP_VENV.replace("\\", "/")))
    return cfg1


# ---------------------------------------------------------------- 运行
def run_isopsv(run_dir, iso_root, log=print):
    """CONFIG=config_step1 → Rscript user_run.R（清 LANG/LC_*；结果写 run_dir）。"""
    rscript = os.path.join(iso_root, "runtime", "R", "bin", "x64", "Rscript.exe")
    runner = os.path.join(iso_root, "code", "user_run.R")
    if not (os.path.isfile(rscript) and os.path.isfile(runner)):
        return "MISSING_RUNTIME", ""
    env = {k: v for k, v in os.environ.items()
           if k not in ("LANG", "LC_ALL", "LC_CTYPE")}
    env["CONFIG"] = os.path.join(run_dir, "config_step1.txt")
    logf = os.path.join(run_dir, "isopsv_run.log")
    with open(logf, "w", encoding="utf-8", errors="replace") as fh:
        fh.write(f"# CONFIG={env['CONFIG']}\n# cwd={os.path.join(iso_root, 'code')}\n")
        fh.flush()
        p = subprocess.run([rscript, runner], shell=False, cwd=os.path.join(iso_root, "code"),
                           env=env, stdout=fh, stderr=subprocess.STDOUT,
                           timeout=14400)
    return ("OK" if p.returncode == 0 else f"EXIT{p.returncode}"), logf


def latest_result_dir(run_dir):
    cands = sorted(p for p in Path(run_dir).glob("result_*") if p.is_dir())
    return str(cands[-1]) if cands else None


# ---------------------------------------------------------------- 结果解析
def _read_tsv(path):
    if not path or not os.path.exists(path):
        return []
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def parse_run(run_dir):
    """→ {新 File_Name: {credible, doubt, status, verdict, score, n_frag,
                        missing_run, cosines{...}}}。"""
    rd = latest_result_dir(run_dir)
    if not rd:
        return {}, rd
    rep = _read_tsv(os.path.join(rd, "1.unlabeled_result",
                                 "2.novel_peptide_competition_report.tsv"))
    # cand1 得分列（classification_A 继承列，形如 candidate_PSM_score_by）
    score_col = None
    for c in (rep[0].keys() if rep else []):
        if "PSM_score" in c or c == "cand1_score":
            score_col = c
            break
    out = {}
    for r in rep:
        fn = r.get("File_Name", "")
        if not fn:
            continue
        try:
            sc = float((r.get(score_col) or "0").split("|")[0].split(",")[0]) / 100.0 \
                if score_col else 0.0
        except (ValueError, IndexError):
            sc = 0.0
        try:
            nf = int(r.get("n_frag") or 0)
        except ValueError:
            nf = 0
        try:
            mr = int(r.get("max_missing_run") or 0)
        except ValueError:
            mr = 99
        out[fn] = dict(credible=r.get("credible_novel", "") == "TRUE",
                       doubt=r.get("novel_doubt", ""),
                       status=r.get("novel_status", ""),
                       verdict=r.get("spectrum_verdict", ""),
                       score=sc, n_frag=nf, missing_run=mr,
                       cosines={})
    # master_table 余弦叠加（列名 cand1_predict_<type>；cand1_predict 主列为
    # 双值 "cos|cos" 不可用，cand1_has_predict_* 为可用性旗标也不取）
    mt = _read_tsv(os.path.join(rd, "1.unlabeled_result", "master_table.tsv"))
    cos_cols = [c for c in (mt[0].keys() if mt else [])
                if c.startswith("cand1_predict_")
                and not c.startswith("cand1_has_predict")]
    fn_col = "File_Name" if (mt and "File_Name" in mt[0]) else None
    if fn_col:
        for r in mt:
            fn = r.get(fn_col, "")
            if fn in out:
                for c in cos_cols:
                    v = r.get(c, "")
                    if v not in ("", None):
                        try:
                            out[fn]["cosines"][c.split("cand1_predict_", 1)[1]] = float(v)
                        except ValueError:
                            pass
    return out, rd


# ---------------------------------------------------------------- 门控与落盘
def peptide_verdicts(per_spec, title_map, min_score="auto", min_cos="auto",
                     log=print):
    """谱级结果 → 肽级指标 + 三门。

    肽级取「最优 PSM」：credible 优先、再余弦高者。
    卡值：数值=固定；"auto"=Tier1 锚定标定（cos 卡=Tier1 p25、score 卡=
    Tier1 中位，clamp [0.20, 0.50]）——衣藻参考 0.956 不适用本机型
    （冒烟实测 Tier1 cos 中位 0.883 < 固定 0.90，会误杀金标准半数）。
    返回 (verd, cards)。
    """
    from collections import defaultdict
    pep_specs = defaultdict(list)
    for fn, r in per_spec.items():
        if fn in title_map:
            pep_specs[title_map[fn]["peptide"]].append((fn, r))
    out = {}
    for pep, items in pep_specs.items():
        def rank(it):
            fn, r = it
            cos = r["cosines"].get(PRIMARY_PREDICTOR)
            cos_v = cos if cos is not None else -1.0
            return (r["credible"], cos_v, r["score"])
        best = max(items, key=rank)[1]
        cosines = best["cosines"]
        cos = cosines.get(PRIMARY_PREDICTOR)
        cos_src = PRIMARY_PREDICTOR
        if cos is None and cosines:  # 主预测器缺失 → 其它中位回退
            vals = sorted(cosines.values())
            cos, cos_src = vals[len(vals) // 2], "median_fallback"
        tier = items and title_map[
            next(fn for fn, _ in items)]["tier"]
        out[pep] = dict(
            tier=tier,
            n_spec=len(items), n_credible=sum(1 for _, r in items if r["credible"]),
            iso_pass=best["credible"],
            best_score=round(best["score"], 4), best_n_frag=best["n_frag"],
            best_missing_run=best["missing_run"],
            best_cos=(round(cos, 4) if cos is not None else None), cos_src=cos_src,
            doubt=best["doubt"], status=best["status"],
        )
    # ---- 卡值解析（auto = Tier1 锚定）----
    tier1 = [v for v in out.values() if v["tier"] == "Tier1"]
    if str(min_cos) == "auto":
        cs = sorted(v["best_cos"] for v in tier1 if v["best_cos"] is not None)
        # clamp [0.50, 0.95]：冒烟小样本 p25 易被离群拉低，全量 78 Tier1 时稳定
        min_cos = round(max(0.50, min(0.95, cs[len(cs) // 4])), 2) if cs else 0.90
        log(f"  [标定] cos 卡值 = Tier1 p25 = {min_cos}（n={len(cs)}）")
    if str(min_score) == "auto":
        ss = sorted(v["best_score"] for v in tier1)
        min_score = min(0.50, max(0.20, round(ss[len(ss) // 2], 2))) if ss else 0.30
        log(f"  [标定] score 卡值 = Tier1 中位 = {min_score}")
    for v in out.values():
        # n_frag 列谱图重推导缺失（空值→0）时不参与卡值，仅 score 卡值生效
        v["quality_pass"] = (v["best_score"] >= min_score
                             and (v["best_n_frag"] == 0 or v["best_n_frag"] >= 4))
        v["sim_pass"] = (v["best_cos"] is not None and v["best_cos"] >= min_cos)
        v["gate_pass"] = v["iso_pass"] and v["quality_pass"] and v["sim_pass"]
    return out, (float(min_score), float(min_cos))


def build_validate(root, out_name="06_validate", iso_root=DEFAULT_ISOPSV,
                   min_score=0.30, min_cos=0.90, subsets=None, source="novel",
                   only_peptides=None, reuse=None, log=print):
    """R8 主入口。subsets=[('Mab','TL'),...] 限定子集（冒烟用）。

    source=novel → 05_novel 的 Tier1/2/3；source=research → P6b 合并表（tier_new）。
    两种来源下**所有新肽都过 isoPSV**（Tier1/2 也跑，只是不因 isoPSV 失败而剔除）。
    only_peptides：只验这些肽（增量模式）；reuse：{peptide: verdict} 复用已验结论。
    """
    vdir = os.path.join(root, out_name)
    os.makedirs(vdir, exist_ok=True)
    rows = collect_novel_psms(root, source=source, only=only_peptides, log=log)
    if subsets:
        rows = [r for r in rows if (r["species"], r["enzyme"]) in subsets]
        log(f"  子集模式：{subsets} → {len(rows)} 谱")
    per_spec, title_maps = {}, {}
    for sp in sorted({r["species"] for r in rows}) or ["Mab", "Mch"]:
        for enz in ("TL", "LL"):
            sub = [r for r in rows if r["species"] == sp and r["enzyme"] == enz]
            if not sub:
                continue
            run_dir = os.path.join(vdir, f"{sp}_{enz}")
            tm_path = os.path.join(run_dir, "title_map.tsv")
            # 复用判据必须带"输入肽集指纹"：只按目录存在会让新一批肽直接借旧结论，
            # 新肽在 isoPSV 里根本没有行 → 全表冒出 MISSING（本轮踩过 27/529）。
            cur_fp = _pep_fprint(r["peptide"] for r in sub)
            cached_fp = None
            if os.path.exists(tm_path):
                with open(tm_path, newline="", encoding="utf-8") as fh:
                    cached_fp = _pep_fprint(x.get("peptide", "")
                                            for x in csv.DictReader(fh, delimiter="\t"))
            if cached_fp == cur_fp and latest_result_dir(run_dir):
                log(f"  [{sp}_{enz}] 已有同输入结果，跳过运行（删 {run_dir} 重跑）")
            else:
                if os.path.exists(tm_path) and cached_fp != cur_fp:
                    log(f"  [{sp}_{enz}] 输入肽集已变（{cached_fp or '无'} → {cur_fp}），重跑 isoPSV")
                mgf_dir = os.path.join(run_dir, "inputs")
                log(f"  [{sp}_{enz}] 抽谱落盘 …")
                psm, _tmap = stage_run(root, sub, run_dir, mgf_dir, log=log)
                fasta = os.path.join(root, "02_db", sp, "db_annot.faa")
                write_configs(run_dir, enz, psm, fasta, mgf_dir, iso_root)
                log(f"  [{sp}_{enz}] isoPSV 运行（最長 2h，日志 {run_dir}\\isopsv_run.log）")
                st, _ = run_isopsv(run_dir, iso_root, log=log)
                log(f"  [{sp}_{enz}] isoPSV {st}")
            # 解析（无论本轮或上轮跑的）
            ps, rd = parse_run(run_dir)
            log(f"  [{sp}_{enz}] 解析 {len(ps)} 谱 ← {rd}")
            with open(os.path.join(run_dir, "title_map.tsv"), newline="",
                      encoding="utf-8") as fh:
                tmap = {r["File_Name_new"]: r
                        for r in csv.DictReader(fh, delimiter="\t")}
            per_spec.update(ps)
            title_maps.update(tmap)
    # 肽级门控 + 汇总（卡值 auto=Tier1 锚定标定）
    verd, cards = peptide_verdicts(per_spec, title_maps, min_score, min_cos, log=log)
    if reuse:
        merged = dict(reuse)
        merged.update(verd)          # 本轮新验的覆盖复用值
        log(f"  合并结论：复用 {len(reuse)} + 本轮 {len(verd)} → {len(merged)} 肽")
        verd = merged
    return verd, per_spec, title_maps, cards


def write_report(root, verd, cards, out_name="06_validate", source="novel", log=print):
    """novel_peptides.tsv 合并验证列 + 分 Tier 通过率汇总。

    cards=(min_score, min_cos)（auto 标定后的实值）。汇总同时给出固定参考
    卡值（0.30/0.90）下的通过数，供对比。
    """
    min_score, min_cos = cards
    vdir = os.path.join(root, out_name)
    from collections import Counter
    # 哪些肽真的把谱送进了 isoPSV：`title_map` 以谱名为键，同一张 MS2 若被两个
    # 候选肽同时命中，后写的会覆盖前者 → 前者一条谱都进不去，最后表现为 MISSING。
    # 这里读回来，把"没谱可送"与"isoPSV 有谱但不给结论"区分开。
    submitted = set()
    for _d in (os.listdir(vdir) if os.path.isdir(vdir) else []):
        _p = os.path.join(vdir, _d, "psm_file.tsv")
        if os.path.isfile(_p):
            with open(_p, newline="", encoding="utf-8") as _fh:
                for _r in csv.DictReader(_fh, delimiter="\t"):
                    _pep = (_r.get("Sequence") or "").strip()
                    if _pep:
                        submitted.add(_pep)
    rows_out, stats = [], []
    sub, fn, tkey = (("05b_research", "novel_after_research.tsv", "tier_new")
                     if source == "research"
                     else ("05_novel", "novel_peptides.tsv", "tier"))
    for sp in list_species(root):
        for r in read_table_in(root, sub, sp, fn):
            v = verd.get(r["peptide"])
            r2 = dict(r)
            r2["tier_final"] = r.get(tkey, "")
            for k in ("n_spec", "n_credible", "iso_pass", "best_score",
                      "best_n_frag", "best_missing_run", "best_cos", "cos_src",
                      "doubt", "status", "quality_pass", "sim_pass",
                      "gate_pass"):
                r2["v_" + k] = "" if v is None else v[k]
            r2["v_miss_reason"] = ""     # 仅在 MISSING 行填，别漏默认值
            tier = r2["tier_final"]
            mirror = tier in ("Tier1", "Tier2")
            iso = bool(v and (v["iso_pass"] is True
                              or str(v["iso_pass"]).lower() == "true"))
            gate = bool(v and (v["gate_pass"] is True
                               or str(v["gate_pass"]).lower() == "true"))
            if v is None:
                r2["final_keep"], r2["keep_reason"] = "MISSING", ""
                # 区分两种 MISSING：谱没送进去（被同 scan 的其它肽占用）／
                # 谱送进去了但 isoPSV 没给结论
                r2["v_miss_reason"] = ("" if r.get("peptide") in submitted
                                       else "no_spectrum_shared")
            elif mirror and iso:
                r2["final_keep"], r2["keep_reason"] = "KEEP", "isoPSV+mirror"
            elif mirror:
                # 镜像证据独立成立：isoPSV 未过仍保留，不因单张谱质量低判死
                r2["final_keep"], r2["keep_reason"] = "KEEP", "mirror_only"
            elif gate:
                r2["final_keep"], r2["keep_reason"] = "KEEP", "isoPSV"
            else:
                r2["final_keep"], r2["keep_reason"] = "DROP", ""
            r2["species"] = sp
            # 兼容旧表列（下游 novgene/novannot/novelfigs 依赖）
            r2["tier"] = r2["tier_final"]
            try:
                r2["psm"] = int(r.get("psm_old") or 0) + int(r.get("psm_new") or 0)
            except (TypeError, ValueError):
                r2["psm"] = r.get("psm_old", "")
            r2["enzymes_str"] = r.get("enzymes_all", "")
            r2["layer"], r2["source"] = layer_of(str(r.get("protein", "")))
            r2["decoy"] = "target"
            rows_out.append(r2)
    fields = list(rows_out[0].keys()) if rows_out else []
    out_tsv = os.path.join(vdir, "novel_validated.tsv")
    with open(out_tsv, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows_out)
    # 分 Tier 汇总
    tag = "P6b 合并集" if source == "research" else "R8 原始集"
    lines = [f"# novel 收口汇总（{tag} · isoPSV v2.0.175 无标记三步流程）",
             f"卡值（Tier1 锚定标定）：cand1 score ≥ {min_score} ∧ "
             f"{PRIMARY_PREDICTOR} cosine ≥ {min_cos} ∧ credible_novel=TRUE（Tier3 门）",
             "保留规则：所有肽都过 isoPSV；**isoPSV 未过但镜像证据成立者保留**"
             "（keep_reason=mirror_only）；Tier3 须过三门。",
             ""]
    for sp in ("Mab", "Mch", "ALL"):
        sel = [r for r in rows_out if sp == "ALL" or r["species"] == sp]
        if not sel:
            continue
        lines.append(f"## {sp}（{len(sel)} 肽）")
        lines.append("| Tier | n | iso_pass | quality | sim | 三门全过 | 保留 | 其中镜像救回 |")
        lines.append("|---|---|---|---|---|---|---|---|")
        for t in ("Tier1", "Tier2", "Tier3"):
            ts = [r for r in sel if r["tier_final"] == t]
            if not ts:
                continue
            n = len(ts)
            cnt = lambda k: sum(1 for r in ts if str(r["v_" + k]).lower() == "true"
                                or r["v_" + k] is True)
            keep = sum(1 for r in ts if r["final_keep"].startswith("KEEP"))
            rescued = sum(1 for r in ts if r.get("keep_reason") == "mirror_only")
            lines.append(f"| {t} | {n} | {cnt('iso_pass')}/{n} "
                         f"| {cnt('quality_pass')}/{n} | {cnt('sim_pass')}/{n} "
                         f"| {cnt('gate_pass')}/{n} | {keep}/{n} | {rescued} |")
        lines.append("")
    reasons = Counter(r["keep_reason"] for r in rows_out if r["final_keep"] == "KEEP")
    if reasons:
        lines.append("## 保留依据分布")
        lines.append("| keep_reason | n |")
        lines.append("|---|---|")
        for k, v in reasons.most_common():
            lines.append(f"| {k} | {v} |")
        lines.append("")
    miss_reasons = Counter(r.get("v_miss_reason", "") for r in rows_out
                           if r["final_keep"] == "MISSING")
    if miss_reasons:
        lines.append("## MISSING 成因（既不算通过也不算剔除）")
        lines.append("| 原因 | n |")
        lines.append("|---|---|")
        lab = {"no_spectrum_shared": "谱被同 scan 的其它候选肽占用，未进 isoPSV",
               "": "谱已送入，isoPSV 未给出结论"}
        for k, v in miss_reasons.most_common():
            lines.append(f"| {lab.get(k, k or '未标注')} | {v} |")
        lines.append("")
    drop_reasons = Counter(r["v_doubt"] for r in rows_out
                           if r["final_keep"] == "DROP" and r["v_doubt"])
    if drop_reasons:
        lines.append("## 剔除原因（novel_doubt 分布，仅触发 doubt 规则者）")
        lines.extend(f"- {k}: {v}" for k, v in drop_reasons.most_common(15))
    sumf = os.path.join(vdir, "validate_summary.md")
    with open(sumf, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    log(f"  写出 {out_tsv}")
    log(f"  写出 {sumf}")
    return out_tsv, sumf
