# -*- coding: utf-8 -*-
"""P6b 候选 ORF 小库回搜（设计见《全流程实验计划》§一·P6b）。

位置：P6 坐标映射之后、P7 镜像匹配/Tier 分级之前。
  1) 取 GSSP 所在 ORF（L2 六框为主，含 L1 Prodigal / L6 转录证据）
  2) 候选 ORF + db_annot（L0+L4+L5）→ 回搜库（远小于 db_6f）
  3) 同一批分 CV MGF 重搜（参数与 E1/E2 完全一致）
  4) 新增肽与新谱图回走同一套四滤（口径不放宽）
  5) 与原始 novel 合并 → 重算 Tier（回搜增量参与分级）

断点续跑：任务目录已有 pFind-Filtered.spectra 即跳过。
子进程安全：参数列表内联 + shell=False（走 pfsearch.run_searcher）。
"""
import csv
import os
import re
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed

from . import pfsearch
from .novel import (digest_krc, digest_krn, iln, isobaric_variants,
                    load_spectra, read_faa)

POOL_FILES = ("L2_6frame.faa", "L1_prodigal_filt.faa", "L1_gms2_filt.faa",
              "L1_augustus_filt.faa", "L6_tx_orfs_filt.faa", "L3_altstart.faa")
# L6 转录 ORF 的名字是 `L6|<contig>|<strand>|<start>-<end>|<n>`（带竖线），
# L3 N 端延伸是 `L3|ext|<protein_id>|up<NN>aa`——都不是下划线前缀。
# 漏了任一条，对应层在 P6b 回搜里会被整层忽略（TP 实测 L3 因此零贡献）。
LAYER_PREFIX = ("L2_", "L1_", "L3_", "L3|", "L6_", "L6|")


# ------------------------------------------------------------------ 小工具
def _core(p):
    """去两端 K/R（镜像核心区）。"""
    while p and p[0] in "KR":
        p = p[1:]
    while p and p[-1] in "KR":
        p = p[:-1]
    return p


def _split_proteins(raw):
    """pFind Proteins 字段可能用 ';' 或 '/' 拼接多蛋白（含 REV_ 诱饵），拆开取层内蛋白。"""
    out = []
    for tok in re.split(r"[;/]", str(raw or "")):
        tok = tok.strip()
        if tok.startswith(LAYER_PREFIX):
            out.append(tok)
    return out


def _orf_locus(protein):
    """L2_<contig>_<frame>_<start> → contig:200bp 分箱；其它层按蛋白 ID。"""
    if protein.startswith("L2_"):
        parts = protein.split("_")
        try:
            start = int(parts[-1])
        except ValueError:
            return protein
        contig = "_".join(parts[1:-2])
        return "%s:%d" % (contig, start // 200 * 200)
    if protein.startswith("L6|"):
        # 与 novel.py 的 tx: locus 口径一致：contig + 200 bp 分箱
        f = protein.split("|")
        contig = f[1] if len(f) > 1 else ""
        try:
            start = int(f[3].split("-")[0])
        except (IndexError, ValueError):
            start = 0
        return "tx:%s:%d" % (contig, start // 200 * 200)
    if protein.startswith("L3"):
        # L3|ext|<protein_id>|up<NN>aa —— 按注释蛋白归并
        f = protein.split("|")
        return "ext:" + (f[2] if len(f) > 2 else protein)
    return "ext:" + protein


def _fw(p):
    return str(p).replace("\\", "/")


# ------------------------------------------------------------------ 任务枚举
def task_index(annot_root, sp):
    """从已有 E1 任务目录反推任务表（保证与 E1/E2 同一批 MGF 与酶）。"""
    base = os.path.join(annot_root, sp)
    out = []
    if not os.path.isdir(base):
        return out
    for short in sorted(os.listdir(base)):
        tdir = os.path.join(base, short)
        if not os.path.isdir(tdir):
            continue
        cfg = os.path.join(tdir, short + ".cfg")
        if not os.path.isfile(cfg):
            continue
        txt = open(cfg, encoding="gbk", errors="replace").read()
        m = re.search(r"^msmspath1=(.+)$", txt, re.M)
        if not m:
            continue
        mgf = m.group(1).strip()
        enz = "LL" if re.search(r"^enzyme=.*LysargiNase", txt, re.M) else "TL"
        out.append(dict(short=short, mgf=mgf, enzyme=enz))
    return out


# ------------------------------------------------------------------ 建库
def candidate_orfs(novel_rows):
    """从 novel 表抽候选 ORF（protein 可能是多蛋白拼接串）。"""
    cand = set()
    for r in novel_rows:
        cand.update(_split_proteins(r.get("protein")))
    return cand


def build_db(root, sp, cand, out_dir, log=print):
    """db_annot + 候选 ORF → 回搜库；返回统计。"""
    db_dir = os.path.join(root, "02_db", sp)
    annot = read_faa(os.path.join(db_dir, "db_annot.faa"))
    have = set(h.split()[0] for h, _ in annot)
    picked, by_pool = [], {}
    for fn in POOL_FILES:
        p = os.path.join(db_dir, fn)
        if not os.path.isfile(p):
            continue
        hit = 0
        for h, s in read_faa(p):
            # 有的源（如 GeneMarkS-2）头里带坐标/注释字段，只取第一个 token 比对，
            # 否则整批 L1_gms2 会被判"未在任何库文件中找到"。
            key = h.split()[0]
            if key in cand and key not in have:
                picked.append((h, s))
                have.add(key)
                hit += 1
        if hit:
            by_pool[fn] = hit
    miss = sorted(cand - set(h.split()[0] for h, _ in picked)
                  - set(h.split()[0] for h, _ in annot))
    os.makedirs(out_dir, exist_ok=True)
    out_faa = os.path.join(out_dir, "db_research_%s.faa" % sp)
    with open(out_faa, "w", encoding="utf-8", newline="\n") as fh:
        for h, s in annot + picked:
            fh.write(">" + h + "\n")
            for i in range(0, len(s), 60):
                fh.write(s[i:i + 60] + "\n")
    n6f = 0
    p6 = os.path.join(db_dir, "db_6f.faa")
    if os.path.isfile(p6):
        n6f = sum(1 for ln in open(p6, encoding="utf-8", errors="replace")
                  if ln.startswith(">"))
    st = dict(species=sp, annot_n=len(annot), cand_n=len(cand),
              picked_n=len(picked), research_n=len(annot) + len(picked),
              db6f_n=n6f, miss_n=len(miss),
              shrink=(n6f / float(len(annot) + len(picked))) if n6f else 0.0,
              by_pool=by_pool)
    log("  [%s] 回搜库：注释 %d + 候选 %d = %d（db_6f %d，缩小 %.1f×）"
        % (sp, st["annot_n"], st["picked_n"], st["research_n"], n6f, st["shrink"]))
    if miss:
        log("  [%s] 警告：%d 个候选 ORF 未在任何库文件中找到（前 3：%s）"
            % (sp, len(miss), "、".join(miss[:3])))
    return out_faa, st


# ------------------------------------------------------------------ 搜库
def _run_one(fasta, tdir_root, t, thread, pfind_bin):
    tdir = os.path.join(tdir_root, t["short"])
    os.makedirs(tdir, exist_ok=True)
    spec = os.path.join(tdir, "pFind-Filtered.spectra")
    if os.path.exists(spec):
        return dict(t, status="reuse", psm=pfsearch.count_result_rows(spec))
    cfg_text = pfsearch.build_cfg_lines(_fw(fasta), _fw(t["mgf"]), _fw(tdir),
                                        t["short"], t["enzyme"], thread=thread)
    cfg = os.path.join(tdir, t["short"] + ".cfg")
    with open(cfg, "w", encoding="gbk", newline="\n") as fh:
        fh.write(cfg_text)
    rc = pfsearch.run_searcher(cfg, pfind_bin)
    n = pfsearch.count_result_rows(spec)
    return dict(t, status="OK" if rc == 0 and n > 0 else "FAIL(rc=%s)" % rc, psm=n)


def run_tasks(fasta, out_root, sp, tasks, pfind_bin, thread=5, workers=4, log=print):
    tr = os.path.join(out_root, sp)
    report, done = [], 0
    with ThreadPoolExecutor(max_workers=max(1, int(workers))) as ex:
        futs = [ex.submit(_run_one, fasta, tr, t, thread, pfind_bin) for t in tasks]
        for f in as_completed(futs):
            r = f.result()
            report.append(r)
            done += 1
            if done % 8 == 0 or done == len(tasks):
                log("  [%s] %d/%d 任务完成" % (sp, done, len(tasks)))
    report.sort(key=lambda x: x["short"])
    return report


# ------------------------------------------------------------------ 收集与判定
def _excluders(root, sp):
    """四滤原料：L0 消化集、L0 串接序列、E1 实测肽。"""
    db_dir = os.path.join(root, "02_db", sp)
    recs = read_faa(os.path.join(db_dir, "db_annot.faa"))
    l0 = [(h, s) for h, s in recs if h.startswith("L0_")]
    known = set()
    for _, s in l0:
        known.update(digest_krc(s))
        known.update(digest_krn(s))
    big = "|".join(s for _, s in l0)
    e1 = set()
    ap = os.path.join(root, "04_search_annot", sp)
    if os.path.isdir(ap):
        for t in os.listdir(ap):
            p = os.path.join(ap, t, "pFind-Filtered.spectra")
            if os.path.exists(p):
                for r in load_spectra(p):
                    e1.add(iln(r.get("Sequence", "")))
    return known, big, e1


def collect(root, sp, out_root, log=print):
    """回搜结果 → 通过同一套四滤的 (peptide, protein) 集合。"""
    known, big, e1 = _excluders(root, sp)
    base = os.path.join(out_root, sp)
    found = {}
    n_raw = n_drop = 0
    for short in sorted(os.listdir(base)) if os.path.isdir(base) else []:
        p = os.path.join(base, short, "pFind-Filtered.spectra")
        if not os.path.isfile(p):
            continue
        enz = "LL" if "_LL_" in short else "TL"
        for r in load_spectra(p):
            prots = _split_proteins(r.get("Proteins"))
            if not prots:
                continue
            prot = prots[0]
            n_raw += 1
            pep = iln(r.get("Sequence", ""))
            if not pep:
                n_drop += 1
                continue
            if pep in known or pep in e1 or pep in big:
                n_drop += 1
                continue
            if any(v in known or v in e1 for v in isobaric_variants(pep) if v != pep):
                n_drop += 1
                continue
            d = found.setdefault((pep, prot),
                                 dict(peptide=pep, protein=prot, psm=0,
                                      enzymes=set(), tasks=set()))
            d["psm"] += 1
            d["enzymes"].add(enz)
            d["tasks"].add(short)
    log("  [%s] 回搜新层命中 PSM %d，过滤后保留肽 %d（剔除 %d）"
        % (sp, n_raw, len(found), n_drop))
    return found


def merge_and_tier(novel_rows, found):
    """原始 novel ∪ 回搜增量 → 重算 Tier。"""
    rows = {}
    for r in novel_rows:
        pep = iln(r["peptide"])
        toks = _split_proteins(r.get("protein")) or [r.get("protein", "")]
        prot = toks[0]
        rows[(pep, prot)] = dict(
            peptide=pep, protein=prot, locus=_orf_locus(prot),
            psm_old=int(r.get("psm") or 0), psm_new=0,
            enz_old=set(x for x in (r.get("enzymes_str") or "").split("/") if x),
            enz_new=set(), tier_old=r.get("tier", ""), is_new=False)
    for (pep, prot), d in found.items():
        if (pep, prot) in rows:
            v = rows[(pep, prot)]
            v["psm_new"] = d["psm"]
            v["enz_new"] = set(d["enzymes"])
        else:
            rows[(pep, prot)] = dict(
                peptide=pep, protein=prot, locus=_orf_locus(prot),
                psm_old=0, psm_new=d["psm"], enz_old=set(),
                enz_new=set(d["enzymes"]), tier_old="", is_new=True)
    core_enz = defaultdict(set)
    for v in rows.values():
        v["enz_all"] = v["enz_old"] | v["enz_new"]
        core_enz[_core(v["peptide"])].update(v["enz_all"])
    t1 = set(c for c, es in core_enz.items() if len(es) >= 2)
    loc_enz = defaultdict(set)
    for v in rows.values():
        loc_enz[v["locus"]].update(v["enz_all"])
    t2 = set(l for l, es in loc_enz.items() if len(es) >= 2)
    for v in rows.values():
        if _core(v["peptide"]) in t1:
            v["tier_new"] = "Tier1"
        elif v["locus"] in t2:
            v["tier_new"] = "Tier2"
        else:
            v["tier_new"] = "Tier3"
        v["enzymes_old"] = "/".join(sorted(v["enz_old"]))
        v["enzymes_new"] = "/".join(sorted(v["enz_new"]))
        v["enzymes_all"] = "/".join(sorted(v["enz_all"]))
    return rows


# ------------------------------------------------------------------ 落盘
def write_species(out_dir, sp, rows, report, dbstat, log=print):
    d = os.path.join(out_dir, sp)
    os.makedirs(d, exist_ok=True)
    # 1) 任务报告
    with open(os.path.join(d, "research_tasks.tsv"), "w", newline="",
              encoding="utf-8-sig") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["task", "enzyme", "status", "psm", "mgf"])
        for r in report:
            w.writerow([r["short"], r["enzyme"], r["status"], r["psm"],
                        os.path.basename(r["mgf"])])
    # 2) 逐肽明细
    fields = ["peptide", "protein", "locus", "tier_old", "tier_new",
              "psm_old", "psm_new", "enzymes_old", "enzymes_new",
              "enzymes_all", "is_new"]
    with open(os.path.join(d, "novel_after_research.tsv"), "w", newline="",
              encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(sorted(rows.values(),
                           key=lambda v: (v["tier_new"], -v["psm_old"],
                                          -v["psm_new"])))
    # 3) Tier 迁移
    mig = Counter()
    for v in rows.values():
        if v["tier_old"] and v["tier_old"] != v["tier_new"]:
            mig[(v["tier_old"], v["tier_new"])] += 1
    with open(os.path.join(d, "tier_migration.tsv"), "w", newline="",
              encoding="utf-8-sig") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["from_tier", "to_tier", "n"])
        for (a, b), n in sorted(mig.items()):
            w.writerow([a, b, n])
    t_new = Counter(v["tier_new"] for v in rows.values())
    t_old = Counter(v["tier_old"] for v in rows.values() if v["tier_old"])
    new_pep = sum(1 for v in rows.values() if v["is_new"])
    psm_gain = sum(v["psm_new"] for v in rows.values())
    psm_old = sum(v["psm_old"] for v in rows.values())
    st = dict(species=sp, db=dbstat, tier_old=dict(t_old), tier_new=dict(t_new),
              n_total=len(rows), n_new_peptide=new_pep,
              psm_old=psm_old, psm_new_research=psm_gain,
              migration=dict(("%s->%s" % k, v) for k, v in mig.items()))
    log("  [%s] 合并后肽 %d（新增 %d）；Tier 旧 %s → 新 %s；新增 PSM %d"
        % (sp, len(rows), new_pep, dict(t_old), dict(t_new), psm_gain))
    return st


def write_summary(out_dir, stats, log=print):
    p = os.path.join(out_dir, "research_summary.md")
    lines = ["# P6b 候选 ORF 小库回搜 · 汇总", "",
             "口径：候选 ORF + db_annot 回搜；同一批分 CV MGF；参数与 E1/E2 一致；",
             "新增肽回走同一套四滤；Tier 由「原始 + 回搜」合并证据重算。", ""]
    lines.append("| 物种 | 注释库 | 候选 ORF | 回搜库 | db_6f | 缩小 | 合并后肽 | 新增肽 | 旧 PSM | 新 PSM |")
    lines.append("|---|---|---|---|---|---|---|---|---|---|")
    for s in stats:
        d = s["db"]
        lines.append("| %s | %d | %d | %d | %d | %.1f× | %d | %d | %d | %d |"
                     % (s["species"], d["annot_n"], d["picked_n"], d["research_n"],
                        d["db6f_n"], d["shrink"], s["n_total"], s["n_new_peptide"],
                        s["psm_old"], s["psm_new_research"]))
    lines += ["", "## Tier 分布（合并证据后）", "",
              "| 物种 | 旧 Tier1/2/3 | 新 Tier1/2/3 | 迁移 |", "|---|---|---|---|"]
    for s in stats:
        o, n = s["tier_old"], s["tier_new"]
        lines.append("| %s | %d/%d/%d | %d/%d/%d | %s |"
                     % (s["species"], o.get("Tier1", 0), o.get("Tier2", 0),
                        o.get("Tier3", 0), n.get("Tier1", 0), n.get("Tier2", 0),
                        n.get("Tier3", 0),
                        ", ".join("%s=%d" % (k, v)
                                  for k, v in sorted(s["migration"].items())) or "无"))
    with open(p, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    log("  写出 %s" % p)
    return p


# ------------------------------------------------------------------ 主入口
def run(workdir, species=("Mab", "Mch"), thread=5, workers=4,
        collect_only=False, log=print):
    root = os.path.abspath(workdir)
    annot_root = os.path.join(root, "04_search_annot")
    out_root = os.path.join(root, "05b_research")
    os.makedirs(out_root, exist_ok=True)
    from . import envdetect
    prow = next((r for r in envdetect.detect_all() if r["name"] == "pFind"), None)
    if not prow or prow.get("status") != "ok":
        raise RuntimeError("pFind 未检测到，无法回搜")
    pfind_bin = os.path.dirname(prow["path"].split(";")[0].strip())
    log("pFind bin: %s" % pfind_bin)
    stats = []
    for sp in species:
        novel_tsv = os.path.join(root, "05_novel", sp, "novel_peptides.tsv")
        if not os.path.isfile(novel_tsv):
            log("  [%s] 缺 %s，跳过" % (sp, novel_tsv))
            continue
        novel_rows = list(csv.DictReader(open(novel_tsv, encoding="utf-8-sig"),
                                         delimiter="\t"))
        cand = candidate_orfs(novel_rows)
        log("  [%s] 原始 novel 肽 %d，涉及候选 ORF %d"
            % (sp, len(novel_rows), len(cand)))
        fasta, dbstat = build_db(root, sp, cand, out_root, log=log)
        tasks = task_index(annot_root, sp)
        log("  [%s] 回搜任务 %d（thread=%d, workers=%d）"
            % (sp, len(tasks), thread, workers))
        if collect_only:
            report = [dict(t, status="skip", psm=pfsearch.count_result_rows(
                os.path.join(out_root, sp, t["short"], "pFind-Filtered.spectra")))
                for t in tasks]
        else:
            report = run_tasks(fasta, out_root, sp, tasks, pfind_bin,
                               thread=thread, workers=workers, log=log)
        found = collect(root, sp, out_root, log=log)
        rows = merge_and_tier(novel_rows, found)
        stats.append(write_species(out_root, sp, rows, report, dbstat, log=log))
    if stats:
        write_summary(out_root, stats, log=log)
    return stats
