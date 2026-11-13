# -*- coding: utf-8 -*-
"""R7 新肽管线：novel 四重过滤 → GSSP（locus 归并）→ Tier 镜像分级。

输入：04_search_annot（E1，注释库结果=排除器②）+ 04_search_disc（E2，发现库结果）
      + 02_db/<sp>/db_annot.faa（L0 蛋白=消化集①与子串滤②原料）
流程（设计书 F 阶段 + 人睾丸线口径）：
  1) E2 中命中 L2_/L3_ 层蛋白的 PSM = novel 候选
  2) 四滤：I/L 归一后 ①∉L0 双酶消化集 ∪ E1 实测靶肽 ②非 L0 蛋白连续≥7aa 子串
     ③等重变体（Q↔E、GG↔N）不落已知集 ④（Trypsin K-P 半特异项不适用于本数据，声明）
  3) 类内 FDR 记账：novel 子集 Target/Decoy 分开计数（pFind 输出列）
  4) GSSP：L2 蛋白按 (contig, 起点 200bp 邻域) 聚 locus；L3 按注释蛋白 ID 归并
  5) Tier（样品级镜像）：Tier1=同一样品组 TL 与 LL 两侧均见（去端 K/R+I/L 归一一致）；
     Tier2=双酶同 locus；Tier3=单酶
"""
import csv
import os
import re
from collections import defaultdict

IL = str.maketrans("IL", "LL")


def iln(s):
    return s.upper().translate(IL)


def digest_krc(seq, min_len=7, max_len=50):
    cuts = [0]
    for i, ch in enumerate(seq[:-1]):
        if ch in "KR" and seq[i + 1] != "P":
            cuts.append(i + 1)
    cuts.append(len(seq))
    return [seq[a:b] for a, b in zip(cuts, cuts[1:]) if min_len <= b - a <= max_len]


def digest_krn(seq, min_len=7, max_len=50):
    cuts = [0]
    for i, ch in enumerate(seq):
        if i > 0 and ch in "KR":
            cuts.append(i)
    cuts.append(len(seq))
    return [seq[a:b] for a, b in zip(cuts, cuts[1:]) if min_len <= b - a <= max_len]


def isobaric_variants(pep):
    """Q↔E、GG↔N 等重变体生成器（有限枚举）。"""
    outs = {pep}
    for i, ch in enumerate(pep):
        if ch == "Q":
            outs.add(pep[:i] + "E" + pep[i + 1:])
        elif ch == "E":
            outs.add(pep[:i] + "Q" + pep[i + 1:])
        if pep[i:i + 2] == "GG":
            outs.add(pep[:i] + "N" + pep[i + 2:])
        if ch == "N":
            for j in range(i + 1, min(i + 3, len(pep))):
                outs.add(pep[:i] + "GG" + pep[j + 1:])
    return outs


def read_faa(path):
    h, buf, out = "", [], []
    for line in open(path, encoding="utf-8", errors="replace"):
        if line.startswith(">"):
            if h:
                out.append((h, "".join(buf)))
            h, buf = line[1:].strip(), []
        else:
            buf.append(line.strip())
    if h:
        out.append((h, "".join(buf)))
    return out


def load_spectra(path):
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def l1_loci(db_dir):
    """L1 ORF → (contig, start1, end1, strand)，来自 `L1_coords.tsv`。

    L1 的蛋白名（`L1_prodigal_6` / `L1_gms2_1`）本身不含坐标，没有这张表就无法
    归 locus，也就算不成同一基因组事件。
    """
    path = os.path.join(db_dir, "L1_coords.tsv")
    out = {}
    if not os.path.exists(path):
        return out
    with open(path, encoding="utf-8-sig") as fh:
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) < 5 or c[0] == "orf":
                continue
            try:
                out[c[0]] = (c[1], int(c[2]), int(c[3]), c[4])
            except ValueError:
                continue
    return out


def build_species_novel(sp, db_dir, annot_root, disc_root, out_dir, log=print):
    """单物种全流程。返回统计 dict。"""
    os.makedirs(out_dir, exist_ok=True)
    l1_map = l1_loci(db_dir)
    # ---- 排除器原料：L0 蛋白与消化集、L0 大序列（子串滤）
    recs = read_faa(os.path.join(db_dir, "db_annot.faa"))
    l0 = [(h, s) for h, s in recs if h.startswith("L0_")]
    known_peps = set()
    for _, s in l0:
        known_peps.update(digest_krc(s))
        known_peps.update(digest_krn(s))
    bigseq = "|".join(s for _, s in l0)  # 分隔符防跨界子串
    log(f"  [{sp}] L0 蛋白 {len(l0)}，消化集 {len(known_peps)} 肽")
    # ---- E1 实测靶肽（排除器②）
    e1_peps = set()
    ap = os.path.join(annot_root, sp)
    for task in os.listdir(ap) if os.path.isdir(ap) else []:
        sp_path = os.path.join(ap, task, "pFind-Filtered.spectra")
        if os.path.exists(sp_path):
            for r in load_spectra(sp_path):
                e1_peps.add(iln(r.get("Sequence", "")))
    log(f"  [{sp}] E1 实测靶肽 {len(e1_peps)}")
    # ---- E2 novel 候选
    cand = []  # (peptide, layer, protein, task, enzyme, decoy)
    dp = os.path.join(disc_root, sp)
    for task in os.listdir(dp) if os.path.isdir(dp) else []:
        sp_path = os.path.join(dp, task, "pFind-Filtered.spectra")
        if not os.path.exists(sp_path):
            continue
        enz = "LL" if "_LL_" in task else "TL"
        for r in load_spectra(sp_path):
            prot = (r.get("Proteins") or "").split(";")[0].rstrip("/")
            # 层白名单：L1 ab initio / L2 六框 / L3 N 端延伸 / L6 转录证据
            # 注意 L3 的蛋白名是 `L3|ext|<protein_id>|up<NN>aa`（竖线），
            # 只认 `L3_` 会把整层静默丢掉（TP 实测漏 30+ 条人工清单里的肽）。
            if not (prot.startswith("L1_") or prot.startswith("L2_")
                    or prot.startswith("L3") or prot.startswith("L6|")):
                continue
            cand.append((r.get("Sequence", ""), prot[:2],
                         prot, task, enz, r.get("Target/Decoy", "")))
    log(f"  [{sp}] E2 层命中 PSM 候选 {len(cand)}")
    # ---- 四滤
    kept = {}
    for seq, layer, prot, task, enz, td in cand:
        p = iln(seq)
        if p in known_peps or p in e1_peps:
            continue                      # ①消化集/实测
        if p in bigseq:
            continue                      # ②≥7aa 连续子串（pFind len≥7 保证）
        if any(v in known_peps or v in e1_peps for v in isobaric_variants(p) if v != p):
            continue                      # ③等重
        key = (p, prot)
        if key not in kept:
            kept[key] = dict(peptide=p, layer=prot[:2], protein=prot, psm=0,
                             enzymes=set(), tasks=set(), decoy=td)
        kept[key]["psm"] += 1
        kept[key]["enzymes"].add(enz)
        kept[key]["tasks"].add(task)
    novs = list(kept.values())
    n_dec = sum(1 for v in novs if v["decoy"].lower().startswith("dec"))
    log(f"  [{sp}] novel 肽 {len(novs)}（含 decoy {n_dec}；类内 decoy 占比 "
        f"{100*n_dec/max(len(novs),1):.2f}%）")
    # ---- GSSP locus：L2 按 (contig, start//200)，L3 按注释蛋白
    locus_of = {}
    for v in novs:
        if v["layer"] == "L1":
            # 与 L2 同口径的 200 bp 分箱（L1 与 L2 命中同一区间即同一事件）
            hit = l1_map.get(v["protein"])
            if hit:
                locus_of[v["peptide"]] = f"{hit[0]}:{hit[1] // 200 * 200}"
            else:
                locus_of[v["peptide"]] = f"ab:{v['protein']}"
            v["source"] = "ab_initio"
        elif v["layer"] == "L2":
            parts = v["protein"].split("_")  # L2_<contig可含_>_frame_start
            start, frame = parts[-1], parts[-2]
            contig = "_".join(parts[1:-2])
            locus_of[v["peptide"]] = f"{contig}:{int(start)//200*200}"
            v["source"] = "intergenic_6frame"
        elif v["layer"] == "L6":
            # L6|<contig>|<strand>|<start>-<end>|_<n>
            f = v["protein"].split("|")
            contig = f[1] if len(f) > 1 else ""
            span = f[3] if len(f) > 3 else ""
            try:
                start = int(span.split("-")[0])
            except ValueError:
                start = 0
            locus_of[v["peptide"]] = f"tx:{contig}:{start//200*200}"
            v["source"] = "transcript_evidence"
        elif v["layer"] == "L3":
            # L3|ext|<protein_id>|up<NN>aa —— 按注释蛋白归 locus（同一基因的不同延伸算一处）
            f = v["protein"].split("|")
            pid = f[2] if len(f) > 2 else v["protein"]
            locus_of[v["peptide"]] = f"ext:{pid}"
            v["source"] = "alt_start_extension"
        else:
            pid = v["protein"].split("_")[2] if len(v["protein"].split("_")) > 2 else v["protein"]
            locus_of[v["peptide"]] = f"ext:{pid}"
            v["source"] = "alt_start_extension"
    # ---- Tier（样品级镜像：TL 与 LL 两侧均见）
    groups = defaultdict(set)  # (species, sample_type) -> 肽（供 Tier1 样品级）
    # 简化实现：peptide 的 enzymes 集合即跨酶证据（同物种内合并）
    def core(p):
        while p and p[0] in "KR":
            p = p[1:]
        while p and p[-1] in "KR":
            p = p[:-1]
        return p
    core_enz = defaultdict(set)
    for v in novs:
        core_enz[core(v["peptide"])].update(v["enzymes"])
    t1_cores = {c for c, es in core_enz.items() if len(es) >= 2}
    t1 = [v for v in novs if core(v["peptide"]) in t1_cores]
    t1_keys = {v["peptide"] for v in t1}
    # Tier2：同 locus 双酶（不同肽）
    loc_enz = defaultdict(set)
    for v in novs:
        loc_enz[locus_of[v["peptide"]]].update(v["enzymes"])
    t2_loci = {loc for loc, es in loc_enz.items() if len(es) >= 2}
    for v in novs:
        if v["peptide"] in t1_keys:
            v["tier"] = "Tier1"
        elif locus_of[v["peptide"]] in t2_loci:
            v["tier"] = "Tier2"
        else:
            v["tier"] = "Tier3"
        v["locus"] = locus_of[v["peptide"]]
        v["enzymes_str"] = "/".join(sorted(v["enzymes"]))
        v["n_tasks"] = len(v["tasks"])
    # ---- 落盘
    fields = ["peptide", "layer", "protein", "source", "locus", "tier",
              "enzymes_str", "psm", "n_tasks", "decoy"]
    with open(os.path.join(out_dir, "novel_peptides.tsv"), "w", newline="",
              encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(sorted(novs, key=lambda x: (x["tier"], -x["psm"])))
    from collections import Counter
    tc = Counter(v["tier"] for v in novs)
    lc = Counter(v["source"] for v in novs)
    stats = dict(species=sp, n_novel=len(novs), n_decoy=n_dec,
                 tier1=tc.get("Tier1", 0), tier2=tc.get("Tier2", 0),
                 tier3=tc.get("Tier3", 0),
                 src_6frame=lc.get("intergenic_6frame", 0),
                 src_altstart=lc.get("alt_start_extension", 0),
                 src_tx=lc.get("transcript_evidence", 0),
                 src_abinitio=lc.get("ab_initio", 0),
                 n_loci=len(loc_enz))
    return stats
