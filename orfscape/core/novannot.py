# -*- coding: utf-8 -*-
"""新编码蛋白：结构预测 / 功能注释 / 同源性比对（深度版）。

结构（序列 + 语言模型两路）：
  · 理化性质   分子量、等电点、GRAVY、不稳定指数、芳香性（Biopython ProtParam）
  · 二级结构   简化 Chou-Fasman（3 残基平滑）+ 置信度
  · 跨膜拓扑   Kyte-Doolittle 19-aa 滑窗 + 正电荷内侧规则 → in/out 拓扑串
  · 卷曲螺旋   七肽重复（a/d 位疏水）滑窗打分
  · 低复杂度   20-aa 窗口最小 Shannon 熵
  · 无序倾向   TOP-IDP 倾向残基占比
  · 接触图     ESM-2（35M）注意力 → 长程接触密度、接触序、接触基 α/β 分数
              （权重缺失时自动跳过并在报告中声明）
功能注释：
  · Swiss-Prot（reviewed）blastp → 策展功能名 / 基因名 / 物种 / 证据等级
  · 本体与交叉物种 RefSeq 蛋白组 blastp → 旁系同源 / 跨种保守性
  · PROSITE 风格基序扫描（糖基化、P-loop、锌指、脂修饰、NLS 等，提示性）
产物：AN_table.tsv/xlsx（annotation + homologs 两表）、AN1–AN4 图
"""
import os
import re
import shutil
import subprocess

import numpy as np

from .. import localpaths
from . import novgene
from .novel import read_faa
from .report import COLORS, COL1, COL2, _apply_zh_font, _save, load_tsv


def _blast_bin():
    """BLAST+ 的 bin 目录：ORFSCAPE_BLAST / paths.local.json 优先，其次 PATH。"""
    p = localpaths.get("blast", r"D:\Bioinformatics\BLAST\ncbi-blast-2.15.0+\bin",
                       r"C:\Program Files\NCBI\blast-2.15.0+\bin")
    if p:
        return p
    hit = shutil.which("makeblastdb") or shutil.which("blastp")
    return os.path.dirname(hit) if hit else ""


BLAST_BIN = _blast_bin()
ANNOT_DIR = os.path.join(os.environ.get("LOCALAPPDATA", os.path.expanduser("~")),
                         "orfscape_annot")
SPROT_GZ = os.path.join(ANNOT_DIR, "uniprot_sprot.fasta.gz")
ESM_PT = os.path.join(ANNOT_DIR, "esm2_t12_35M.pt")

KD = dict(A=1.8, R=-4.5, N=-3.5, D=-3.5, C=2.5, Q=-3.5, E=-3.5, G=-0.4,
          H=-3.2, I=4.5, L=3.8, K=-3.9, M=1.9, F=2.8, P=-1.6, S=-0.8,
          T=-0.7, W=-0.9, Y=-1.3, V=4.2)
PA = dict(A=1.42, R=0.98, N=0.67, D=1.01, C=0.70, Q=1.11, E=1.51, G=0.57,
          H=1.00, I=1.08, L=1.21, K=1.16, M=1.45, F=1.13, P=0.57, S=0.77,
          T=0.83, W=1.08, Y=0.69, V=1.06)
PB = dict(A=0.83, R=0.93, N=0.89, D=0.54, C=1.19, Q=1.10, E=0.37, G=0.75,
          H=0.87, I=1.60, L=1.30, K=0.74, M=1.05, F=1.38, P=0.55, S=0.75,
          T=1.19, W=1.37, Y=1.47, V=1.70)
DISORDER = set("AEKQRSTPGHNDC")
MOTIFS = [
    ("N-glycosylation", r"N[^P][ST][^P]"),
    ("P-loop (NTP bind)", r"[AG][AGS][A-Z]{2}GKS"),
    ("Zinc finger C2H2", r"C[A-Z]{2,4}C[A-Z]{3}[LIVMFY][A-Z]{8}H[A-Z]{3,5}H"),
    ("N-myristoylation", r"^G[^EDRKHPFYW][STAGCN][^P]"),
    ("CaaX prenylation", r"C[ACD][ACDILMVF][ACDGKNRS]$"),
    ("Bipartite NLS", r"^[KR][KR][A-Z]{10,12}[KR]{3}"),
    ("EF-hand", r"D[DNQSTAP][DNQGHRK]G[LIVMC][DENQSTAGC][A-Z]{2}[DE]"),
    ("NES-like", r"L[A-Z]{2,3}[LIVFM][A-Z]{2,3}L[A-Z]L"),
    ("His-Asp-Ser triad", r"H[DE][A-Z]{2}[ST]"),
    ("ATP/GTP site A", r"[AG]x{4}GKS"),
]


def T(zh, en):
    from .i18n import get_lang
    return zh if get_lang() == "zh" else en


# ---------------------------------------------------------------- 目标蛋白
def novel_orfs(root):
    """新编码蛋白（基因间区类、可信新肽所属 ORF）→ ({orf: seq}, {orf: species})。"""
    vdir = os.environ.get("ORFSCAPE_VALIDATE_DIR", "06_validate")
    src = os.path.join(root, vdir, "novel_validated.tsv")
    rows = load_tsv(src) if os.path.exists(src) else []
    rows = [r for r in rows
            if str(r.get("final_keep", "")).startswith("KEEP")]
    rows = novgene.classify_events(root, rows, log=lambda *a: None)
    keep = {r["event_orf"] for r in rows
            if r.get("event_type") == novgene.TYPE_INTERGENIC}
    seqs, sp_of = {}, {}
    for sp in novgene._species_list(rows):
        # L2 六框 + L6 转录 ORF 两处都可能承载新编码基因
        for fn in ("L2_6frame.faa", "L6_tx_orfs_filt.faa"):
            faa = os.path.join(root, "02_db", sp, fn)
            if not os.path.exists(faa):
                continue
            for h, s in read_faa(faa):
                name = h.split()[0]
                if name in keep:
                    seqs[name] = s
                    sp_of[name] = sp
    return seqs, sp_of


# ---------------------------------------------------------------- 结构特征
def _clean_seq(seq):
    """去掉终止符与模糊字符：Prodigal 的 -a 输出会在 C 端带 '*'，
    Biopython ProtParam 遇到 '*/X/U' 会直接抛 ValueError。"""
    s = str(seq or "").rstrip("*").upper()
    return "".join(ch for ch in s if ch in "ACDEFGHIKLMNPQRSTVWY")


def physchem(seq):
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
    a = ProteinAnalysis(_clean_seq(seq))
    return dict(mw=a.molecular_weight(), pI=a.isoelectric_point(),
                gravy=a.gravy(), instability=a.instability_index(),
                aromaticity=a.aromaticity())


def ss_chou_fasman(seq):
    al = np.convolve([PA.get(a, 1.0) for a in seq], np.ones(3) / 3, "same")
    be = np.convolve([PB.get(a, 1.0) for a in seq], np.ones(3) / 3, "same")
    helix = int(((al >= be) & (al > 1.0)).sum())
    beta = int(((be > al) & (be > 1.0)).sum())
    n = max(len(seq), 1)
    frac = (helix / n, beta / n, (n - helix - beta) / n)
    conf = float(np.mean(np.abs(al - be))) if len(seq) else 0.0
    if frac[0] >= 0.5:
        cls = "all-alpha"
    elif frac[1] >= 0.5:
        cls = "all-beta"
    elif frac[0] >= 0.25 and frac[1] >= 0.25:
        cls = "alpha/beta"
    else:
        cls = "mixed/coil-rich"
    return frac, cls, conf


def tm_topology(seq, win=19, thr=1.6):
    """KD 滑窗找跨膜段 + 正电荷内侧规则 → (段数, 段列表, 拓扑串)。"""
    if len(seq) < win:
        return 0, [], ""
    v = np.array([KD.get(a, 0.0) for a in seq])
    raw, in_seg, st = [], False, 0
    for i in range(len(v) - win + 1):
        m = v[i:i + win].mean()
        if m >= thr and not in_seg:
            in_seg, st = True, i
        elif m < thr and in_seg:
            in_seg = False
            raw.append((st + 1, i + win - 1))
    if in_seg:
        raw.append((st + 1, len(seq)))
    segs = []
    for s, e in raw:
        if segs and s - segs[-1][1] <= 3:
            segs[-1] = (segs[-1][0], e)
        else:
            segs.append((s, e))
    topo, last = [], 0
    for i, (s, e) in enumerate(segs):
        loop = seq[last:s - 1]
        pos = sum(1 for a in loop if a in "KR")
        neg = sum(1 for a in loop if a in "DE")
        side = "in" if i % 2 == 0 else "out"
        if i > 0 and pos - neg < 0:
            side = "out" if side == "in" else "in"
        topo.append(side)
        last = e
    return len(segs), segs, "-".join(topo)


def coiled_coil(seq, win=28, thr=0.5):
    hyd = set("AVILMFWY")
    hits = tot = 0
    for i in range(max(len(seq) - win + 1, 1)):
        w = seq[i:i + win]
        if len(w) < win:
            break
        tot += 1
        score = 0.0
        for k in (0, 3):
            pos = [w[j] for j in range(k, len(w), 7)]
            if pos:
                score += sum(1 for a in pos if a in hyd) / len(pos)
        if score / 2 >= thr:
            hits += 1
    return hits / max(tot, 1)


def low_complexity(seq, win=20):
    if len(seq) < win:
        return float("nan")
    ents = []
    for i in range(len(seq) - win + 1):
        w = seq[i:i + win]
        f = np.array([w.count(a) for a in set(w)], dtype=float) / len(w)
        ents.append(float(-(f * np.log2(f)).sum()))
    return float(min(ents))


def disorder_frac(seq):
    return sum(1 for a in seq if a in DISORDER) / max(len(seq), 1)


def motif_scan(seq):
    return [name for name, pat in MOTIFS if re.search(pat, seq)]


# ---------------------------------------------------------------- ESM 接触图
_ESM = {"tried": False, "model": None}


def _esm(log=print):
    if _ESM["tried"]:
        return _ESM["model"]
    _ESM["tried"] = True
    try:
        import esm
        import esm.pretrained as P
        # 本地只有主权重（无 contact-regression 头）：跳过其加载，
        # 仅用注意力图作粗结构信号（报告中已声明未做监督校准）
        P._has_regression_weights = lambda name: False
        m, alphabet = P.load_model_and_alphabet_local(ESM_PT)
        m.eval()
        _ESM["model"] = (m, alphabet)
    except Exception as e:
        log(f"  [ESM] 权重不可用（{type(e).__name__}），接触图跳过")
        _ESM["model"] = None
    return _ESM["model"]


def _esm_cache_path(orf):
    d = os.path.join(ANNOT_DIR, "esm_cache")
    os.makedirs(d, exist_ok=True)
    return os.path.join(d, orf.replace("/", "_") + ".npz")


def esm_contacts(seq, log=print, orf=None):
    """ESM-2 注意力接触图统计；按 ORF 名缓存（npz）避免重算。"""
    cpath = _esm_cache_path(orf) if orf else None
    if cpath and os.path.exists(cpath):
        try:
            z = np.load(cpath)
            return dict(lr_contact=float(z["lr"]),
                        contact_order=float(z["co"]),
                        esm_helix=float(z["helix"]),
                        esm_beta=float(z["beta"]), A=z["A"])
        except Exception:
            pass
    got = _esm(log=log)
    if got is None:
        return None
    m, alphabet = got
    try:
        import torch
        _, _, toks = alphabet.get_batch_converter()([("p", seq)])
        with torch.no_grad():
            out = m(toks, repr_layers=[], need_head_weights=True)
            att = out["attentions"]
            am = att[-1].mean(dim=1)[0].numpy()
        L = len(seq)
        A = am[1:L + 1, 1:L + 1]
        A = 0.5 * (A + A.T)
        k = max(int(round(L / 5)), 5)
        C = np.zeros((L, L), dtype=bool)
        for i in range(L):
            order = np.argsort(-A[i])
            picked = [int(j) for j in order if abs(int(j) - i) >= 6][:k]
            for j in picked:
                C[i, j] = C[j, i] = True
        lr = float(C.sum()) / max(L * L, 1)
        ii, jj = np.where(C)
        co = float(np.mean(np.abs(ii - jj)) / max(L, 1)) if len(ii) else 0.0
        d = np.abs(np.subtract.outer(np.arange(L), np.arange(L)))
        total = max(int(C.sum()), 1)
        helix = float((C & ((d == 3) | (d == 4))).sum()) / total
        beta = float((C & (d == 2)).sum()) / total
        if cpath:
            try:
                np.savez_compressed(cpath, lr=lr, co=co, helix=helix,
                                    beta=beta, A=A.astype(np.float32))
            except Exception:
                pass
        return dict(lr_contact=lr, contact_order=co, esm_helix=helix,
                    esm_beta=beta, A=A)
    except Exception as e:
        log(f"  [ESM] 单蛋白失败：{type(e).__name__}")
        return None


# ---------------------------------------------------------------- 同源
def _make_db(faa, tag, log=print):
    cache = os.path.join(ANNOT_DIR, "blastdb")
    os.makedirs(cache, exist_ok=True)
    base = os.path.join(cache, "db_" + tag)
    if os.path.exists(base + ".pdb"):
        return base
    mk = os.path.join(BLAST_BIN, "makeblastdb.exe")
    if not os.path.exists(mk):
        return None
    r = subprocess.run([mk, "-in", faa, "-dbtype", "prot", "-out", base],
                       shell=False, capture_output=True)
    if not os.path.exists(base + ".pdb"):
        if r.stderr:
            log("  [blast] 建库失败：" + r.stderr.decode("utf-8", "replace")[-160:])
        return None
    return base


def _blast(query, db, out_path, log=print, top=5):
    bp = os.path.join(BLAST_BIN, "blastp.exe")
    if not db or not os.path.exists(bp):
        return []
    cmd = [bp, "-query", query, "-db", db, "-evalue", "1e-3",
           "-max_target_seqs", str(top), "-num_threads", "8",
           "-outfmt", "6 qseqid sseqid pident length qlen slen evalue stitle",
           "-out", out_path]
    r = subprocess.run(cmd, shell=False, capture_output=True)
    if r.stderr:
        msg = r.stderr.decode("utf-8", "replace")
        if "Error" in msg:
            log("  [blast] " + msg.strip()[-160:])
    hits = []
    if os.path.exists(out_path):
        with open(out_path, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 8:
                    continue
                hits.append(dict(q=f[0], s=f[1], pid=float(f[2]),
                                 aln=int(f[3]), qlen=int(f[4]),
                                 slen=int(f[5]), ev=float(f[6]), title=f[7]))
    return hits


def sprot_meta(title):
    t = title or ""
    name = t
    m = re.match(r"^sp\|[^|]+\|[^ ]+\s+(.*)$", t)
    if m:
        name = m.group(1)
    org = ""
    mo = re.search(r"\sOS=([^=]+?)(?:\s+OX=|\s+GN=|\s+PE=|\s+SV=|$)", t)
    if mo:
        org = mo.group(1).strip()
    gene = ""
    mg = re.search(r"\sGN=([^\s]+)", t)
    if mg:
        gene = mg.group(1)
    ev = ""
    me = re.search(r"\sPE=(\d)", t)
    if me:
        ev = me.group(1)
    name = re.sub(r"\s(OS|OX|GN|PE|SV)=.*$", "", name).strip()
    return name[:70], gene, org[:50], ev


def homology_class(pid, cov):
    if pid is None:
        return "none"
    if pid >= 99 and cov >= 99:
        return "identical"
    if pid >= 70 and cov >= 70:
        return "strong"
    if cov >= 50 and pid >= 30:
        return "moderate"
    return "weak"


# ---------------------------------------------------------------- 汇总
def annotate_all(seqs, sp_of, use_esm=True, log=print):
    ann = []
    for i, (name, seq) in enumerate(sorted(seqs.items()), 1):
        pc = physchem(seq)
        ss, cls, conf = ss_chou_fasman(seq)
        tm_n, tm_segs, topo = tm_topology(seq)
        rec = dict(orf=name, species=sp_of.get(name, ""), length=len(seq),
                   mw=pc["mw"], pI=pc["pI"], gravy=pc["gravy"],
                   instability=pc["instability"], aromaticity=pc["aromaticity"],
                   ss=ss, ss_cls=cls, ss_conf=conf, tm=tm_n,
                   tm_segs=tm_segs, topology=topo,
                   coiled=coiled_coil(seq), lc=low_complexity(seq),
                   disorder=disorder_frac(seq), motifs=motif_scan(seq),
                   homology="none")
        if use_esm:
            e = esm_contacts(seq, log=log, orf=name)
            if e:
                rec.update(lr_contact=e["lr_contact"],
                           contact_order=e["contact_order"],
                           esm_helix=e["esm_helix"], esm_beta=e["esm_beta"],
                           _A=e["A"])
        ann.append(rec)
        if i % 25 == 0:
            log(f"    结构预测 {i}/{len(seqs)}")
    return ann


def run_homology(root, ann, seqs, out_dir, log=print):
    qdir = os.path.join(out_dir, "blast")
    os.makedirs(qdir, exist_ok=True)
    by_sp = {}
    for r in ann:
        by_sp.setdefault(r["species"], []).append(r)
    db_sp = None
    if os.path.exists(SPROT_GZ):
        fasta = os.path.join(ANNOT_DIR, "uniprot_sprot.fasta")
        if not os.path.exists(fasta):
            import gzip
            with gzip.open(SPROT_GZ, "rb") as fi, open(fasta, "wb") as fo:
                for chunk in iter(lambda: fi.read(1 << 20), b""):
                    fo.write(chunk)
        db_sp = _make_db(fasta, "sprot", log=log)
    for sp, items in by_sp.items():
        qfaa = os.path.join(qdir, f"novel_{sp}.faa")
        with open(qfaa, "w", encoding="utf-8") as fh:
            for r in items:
                fh.write(">" + r["orf"] + chr(10))
                s = seqs[r["orf"]]
                for i in range(0, len(s), 60):
                    fh.write(s[i:i + 60] + chr(10))
        hits = (_blast(qfaa, db_sp, os.path.join(qdir, f"sprot_{sp}.tsv"),
                       log=log) if db_sp else [])
        best = {}
        for h in hits:
            cov = 100.0 * h["aln"] / max(h["qlen"], 1)
            if h["q"] not in best or h["pid"] * cov > best[h["q"]]["pid"] * best[h["q"]]["cov"]:
                h["cov"] = cov
                best[h["q"]] = h
        for r in items:
            b = best.get(r["orf"])
            if b:
                nm, gene, org, ev = sprot_meta(b["title"])
                r.update(sp_hit=b["s"], sp_ident=b["pid"], sp_cov=b["cov"],
                         sp_eval=b["ev"], sp_name=nm, sp_gene=gene,
                         sp_org=org, sp_ev=ev)
            else:
                r.update(sp_hit="", sp_ident=None, sp_cov=None, sp_eval=None,
                         sp_name="", sp_gene="", sp_org="", sp_ev="")
        # 跨物种 RefSeq 对照：只比"确实有基因组目录"的物种（单物种项目就没有 cross）
        _others = [x for x in os.listdir(os.path.join(root, "01_genome"))
                   if os.path.isdir(os.path.join(root, "01_genome", x))]
        for other in ([sp] + [x for x in sorted(_others) if x != sp]):
            faa = os.path.join(root, "01_genome", other, "protein.faa")
            if not os.path.exists(faa):
                continue
            db2 = _make_db(faa, other, log=log)
            hits2 = _blast(qfaa, db2, os.path.join(qdir, f"ref_{other}_{sp}.tsv"),
                           log=log)
            best2 = {}
            for h in hits2:
                cov = 100.0 * h["aln"] / max(h["qlen"], 1)
                if h["q"] not in best2 or h["pid"] * cov > best2[h["q"]]["pid"] * best2[h["q"]]["cov"]:
                    h["cov"] = cov
                    best2[h["q"]] = h
            key = "self" if other == sp else "cross"
            for r in items:
                b = best2.get(r["orf"])
                if b:
                    fn = re.sub(r"^[A-Z]{2}_\d+\.\d+\s*", "", b["title"])
                    r.update({f"{key}_hit": b["s"], f"{key}_ident": b["pid"],
                              f"{key}_cov": b["cov"], f"{key}_func": fn[:70]})
                else:
                    r.update({f"{key}_hit": "", f"{key}_ident": None,
                              f"{key}_cov": None, f"{key}_func": ""})
    for r in ann:
        pid, cov, src = r.get("sp_ident"), r.get("sp_cov"), "Swiss-Prot"
        if pid is None:
            pid, cov = r.get("self_ident"), r.get("self_cov")
            src = "RefSeq-self" if pid is not None else ""
        r["homology"] = homology_class(pid, cov)
        r["homology_src"] = src
        r["function"] = r.get("sp_name") or r.get("self_func") or ""
    return ann


# ---------------------------------------------------------------- 输出
def _col_letter(n):
    """1→A, 26→Z, 27→AA。"""
    s = ""
    while n > 0:
        n, r = divmod(n - 1, 26)
        s = chr(65 + r) + s
    return s


def write_tables(out_dir, ann, zh):
    def f(v, nd=2):
        if v is None or v == "":
            return ""
        return f"{v:.{nd}f}" if isinstance(v, float) else str(v)

    heads = [T("ORF", "ORF"), T("物种", "Species"), T("长度", "Length"),
             T("MW(kDa)", "MW(kDa)"), T("pI", "pI"), T("GRAVY", "GRAVY"),
             T("不稳定", "Instab."), T("α%", "α%"), T("β%", "β%"),
             T("SS 类别", "SS class"), T("SS 置信", "SS conf"),
             T("TM 段", "TM"), T("拓扑", "Topology"),
             T("卷曲螺旋", "Coiled"), T("低复杂熵", "LC entropy"),
             T("无序", "Disorder"), T("长程接触", "LR contact"),
             T("接触序", "Contact order"), T("ESM α", "ESM a"),
             T("ESM β", "ESM b"), T("基序", "Motifs"),
             T("同源等级", "Homology"), T("同源来源", "Source"),
             T("SP 命中", "SP hit"), T("SP 一致性", "SP ident"),
             T("SP 覆盖", "SP cov"), T("SP 蛋白名", "SP protein"),
             T("SP 基因", "SP gene"), T("SP 物种", "SP organism"),
             T("SP 证据", "SP PE"), T("RefSeq 本体", "RefSeq self"),
             T("交叉物种", "Cross-species")]
    rows = [heads]
    for r in ann:
        rows.append([
            r["orf"], r["species"], str(r["length"]), f(r["mw"] / 1000, 1),
            f(r["pI"]), f(r["gravy"]), f(r["instability"], 1),
            f"{100*r['ss'][0]:.0f}", f"{100*r['ss'][1]:.0f}", r["ss_cls"],
            f(r.get("ss_conf")), str(r["tm"]), r.get("topology", ""),
            f(r.get("coiled")), f(r.get("lc")), f(r["disorder"]),
            f(r.get("lr_contact"), 4), f(r.get("contact_order"), 3),
            f(r.get("esm_helix"), 3), f(r.get("esm_beta"), 3),
            "; ".join(r.get("motifs", [])), r["homology"],
            r.get("homology_src", ""), r.get("sp_hit", ""),
            f(r.get("sp_ident"), 1), f(r.get("sp_cov"), 0),
            r.get("sp_name", ""), r.get("sp_gene", ""), r.get("sp_org", ""),
            r.get("sp_ev", ""), r.get("self_func", ""), r.get("cross_func", ""),
        ])
    tbl = os.path.join(out_dir, "AN_table.tsv")
    with open(tbl, "w", newline="", encoding="utf-8-sig") as fh:
        for r in rows:
            fh.write("\t".join(r) + "\n")
    try:
        from openpyxl import Workbook
        from openpyxl.styles import Alignment, Font
        wb = Workbook()
        ws = wb.active
        ws.title = "annotation"
        for r in rows:
            ws.append(r)
        for c in ws[1]:
            c.font = Font(name="Arial", size=9, bold=True)
            c.alignment = Alignment(wrap_text=True, vertical="top")
        for row in ws.iter_rows(min_row=2):
            for c in row:
                c.font = Font(name="Arial", size=9)
        widths = [22, 7, 7, 9, 6, 7, 7, 6, 6, 13, 7, 6, 10, 8, 9, 7, 9, 9,
                  6, 6, 26, 10, 11, 14, 9, 8, 26, 12, 20, 6, 26, 26]
        for i, w in enumerate(widths, 1):
            ws.column_dimensions[_col_letter(i)].width = w
        ws2 = wb.create_sheet("homologs")
        ws2.append(heads)
        for c in ws2[1]:
            c.font = Font(name="Arial", size=9, bold=True)
        for r in rows[1:]:
            if r[21] != "none":
                ws2.append(r)
        book = os.path.join(out_dir, "AN_table.xlsx")
        wb.save(book)
    except ImportError:
        pass
    return tbl


def plot_properties(out_dir, ann, zh):
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 3, figsize=(COL2, 4.4),
                             gridspec_kw={"wspace": 0.32, "hspace": 0.60})
    specs = [("length", T("长度 (aa)", "Length (aa)"), 1,
              T("a  长度", "a  Length")),
             ("mw", T("分子量 (kDa)", "MW (kDa)"), 1e-3, T("b  分子量", "b  MW")),
             ("pI", "pI", 1, T("c  等电点", "c  pI")),
             ("gravy", "GRAVY", 1, T("d  疏水性", "d  GRAVY")),
             ("disorder", T("无序比例", "Disorder frac"), 1,
              T("e  无序倾向", "e  Disorder")),
             ("lr_contact", T("长程接触密度", "LR contact density"), 1,
              T("f  ESM 接触密度", "f  ESM contacts"))]
    for ax, (key, xlab, scale, ttl) in zip(axes.ravel(), specs):
        vals = [r[key] * scale for r in ann if r.get(key) is not None]
        if not vals:
            ax.axis("off")
            ax.set_title(ttl, fontsize=7, fontweight="bold", loc="left")
            ax.text(0.5, 0.5, T("无数据", "no data"), ha="center",
                    transform=ax.transAxes, fontsize=5, color="#9CA3AF")
            continue
        ax.hist(vals, bins=20, color=COLORS[0], alpha=0.85)
        med = float(np.median(vals))
        ax.axvline(med, color=COLORS[3], ls="--", lw=0.8)
        ax.set_xlabel(xlab, fontsize=5.6)
        ax.set_ylabel(T("蛋白数", "Proteins"), fontsize=5.6)
        ax.set_title(ttl, fontsize=7, fontweight="bold", loc="left")
        ax.text(med, ax.get_ylim()[1] * 0.95, f" med {med:.2f}",
                fontsize=4.6, color=COLORS[3], va="top")
    _save(fig, out_dir, "AN1_properties")


def plot_homology(out_dir, ann, zh):
    import matplotlib.pyplot as plt
    from collections import Counter
    pts = [(r.get("sp_ident"), r.get("sp_cov")) for r in ann
           if r.get("sp_ident") is not None]
    fig, axes = plt.subplots(1, 2, figsize=(COL2, 2.8),
                             gridspec_kw={"wspace": 0.34})
    ax = axes[0]
    if pts:
        ax.scatter([p[0] for p in pts], [p[1] for p in pts], s=14,
                   color=COLORS[0], alpha=0.8, edgecolors="white", linewidths=0.3)
        ax.axhline(70, color="#9CA3AF", ls=":", lw=0.7)
        ax.axvline(70, color="#9CA3AF", ls=":", lw=0.7)
    else:
        ax.text(0.5, 0.5, T("Swiss-Prot 库未就绪", "Swiss-Prot DB not ready"),
                ha="center", transform=ax.transAxes, fontsize=5.5,
                color="#9CA3AF")
    ax.set_xlabel(T("一致性 (%)", "Identity (%)"), fontsize=5.6)
    ax.set_ylabel(T("覆盖度 (%)", "Coverage (%)"), fontsize=5.6)
    ax.set_title(T("a  Swiss-Prot 命中", "a  Swiss-Prot hits"), fontsize=7,
                 fontweight="bold", loc="left")
    ax = axes[1]
    order = ["identical", "strong", "moderate", "weak", "none"]
    labels = {"identical": T("一致", "Identical"), "strong": T("强", "Strong"),
              "moderate": T("中等", "Moderate"), "weak": T("弱", "Weak"),
              "none": T("无", "None")}
    cnt = Counter(r["homology"] for r in ann)
    vals = [cnt.get(k, 0) for k in order]
    bars = ax.bar(range(len(order)), vals, 0.6,
                  color=[COLORS[3], COLORS[1], COLORS[2], COLORS[0], "#9CA3AF"])
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2, v + max(vals + [1]) * 0.02,
                f"{v}", ha="center", fontsize=5)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([labels[k] for k in order], fontsize=5, rotation=20,
                       ha="right")
    ax.set_ylabel(T("蛋白数", "Proteins"), fontsize=5.6)
    ax.set_title(T("b  同源等级", "b  Homology classes"), fontsize=7,
                 fontweight="bold", loc="left")
    _save(fig, out_dir, "AN2_homology")


def plot_features(out_dir, ann, zh):
    import matplotlib.pyplot as plt
    from collections import Counter
    fig, axes = plt.subplots(1, 3, figsize=(COL2, 2.4),
                             gridspec_kw={"wspace": 0.38})
    ax = axes[0]
    tm = Counter(r["tm"] for r in ann)
    ks = sorted(tm)
    ax.bar([str(k) for k in ks], [tm[k] for k in ks], 0.6, color=COLORS[0])
    ax.set_xlabel(T("跨膜段数", "TM segments"), fontsize=5.6)
    ax.set_ylabel(T("蛋白数", "Proteins"), fontsize=5.6)
    ax.set_title(T("a  跨膜螺旋", "a  TM helices"), fontsize=7,
                 fontweight="bold", loc="left")
    ax = axes[1]
    ss = Counter(r["ss_cls"] for r in ann)
    ks = list(ss)
    ax.bar(range(len(ks)), [ss[k] for k in ks], 0.6, color=COLORS[1])
    ax.set_xticks(range(len(ks)))
    ax.set_xticklabels(ks, fontsize=5, rotation=20, ha="right")
    ax.set_ylabel(T("蛋白数", "Proteins"), fontsize=5.6)
    ax.set_title(T("b  二级结构类别", "b  SS class"), fontsize=7,
                 fontweight="bold", loc="left")
    ax = axes[2]
    mf = Counter(m for r in ann for m in r.get("motifs", []))
    if mf:
        items = mf.most_common(8)
        ax.barh([k for k, _ in items][::-1], [v for _, v in items][::-1],
                0.6, color=COLORS[2])
        ax.set_xlabel(T("蛋白数", "Proteins"), fontsize=5.6)
        ax.tick_params(labelsize=4.8)
    ax.set_title(T("c  基序扫描（提示性）", "c  Motif hits"), fontsize=7,
                 fontweight="bold", loc="left")
    _save(fig, out_dir, "AN3_features")


def plot_contacts(out_dir, ann, zh, top=6):
    import matplotlib.pyplot as plt
    cand = [r for r in ann if r.get("_A") is not None]
    if not cand:
        return
    cand = sorted(cand, key=lambda r: (r["homology"] == "none", r["length"]),
                  reverse=True)[:top]
    cols = 3
    rows_n = int(np.ceil(len(cand) / cols))
    fig, axes = plt.subplots(rows_n, cols, figsize=(COL2, 1.9 * rows_n + 0.4),
                             gridspec_kw={"wspace": 0.28, "hspace": 0.5})
    axes = np.atleast_1d(axes).ravel()
    for ax, r in zip(axes, cand):
        ax.imshow(r["_A"], cmap="viridis", interpolation="nearest")
        ax.set_title(f"{r['orf'][-14:]} ({r['length']} aa)", fontsize=4.8)
        ax.tick_params(labelsize=3.6)
    for ax in axes[len(cand):]:
        ax.axis("off")
    fig.suptitle(T("ESM-2 注意力接触图（示例）",
                   "ESM-2 attention contact maps (examples)"),
                 fontsize=7.4, fontweight="bold", x=0.01, ha="left", y=1.02)
    _save(fig, out_dir, "AN4_contact_maps")


# ---------------------------------------------------------------- 入口
def build_novannot(root, out_dir, use_esm=True, log=print):
    from .i18n import get_lang
    zh = get_lang() == "zh"
    _apply_zh_font()
    os.makedirs(out_dir, exist_ok=True)
    seqs, sp_of = novel_orfs(root)
    log(f"  新编码蛋白 {len(seqs)} 个 → 结构预测（ESM={'on' if use_esm else 'off'}）")
    ann = annotate_all(seqs, sp_of, use_esm=use_esm, log=log)
    log("  同源比对（Swiss-Prot + RefSeq 本体/交叉）")
    ann = run_homology(root, ann, seqs, out_dir, log=log)
    write_tables(out_dir, ann, zh)
    plot_properties(out_dir, ann, zh)
    plot_homology(out_dir, ann, zh)
    plot_features(out_dir, ann, zh)
    if use_esm:
        plot_contacts(out_dir, ann, zh)
    from collections import Counter
    log(f"  完成：同源 {dict(Counter(r['homology'] for r in ann))}")
    return ann
