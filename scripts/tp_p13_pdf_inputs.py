# -*- coding: utf-8 -*-
"""TP P13 第 4 件：为新版 PDF 引擎准备输入（Table1 / Table3 / spectra / cache / predict）。

把 ORFscape 的 ORF 命名解析成基因组坐标：
  L2_<contig>_<frame><off>_<start>   → 六框 ORF（长度取 L2_6frame.faa）
  L6|<contig>|<strand>|<s>-<e>|<n>   → 转录 ORF（名字自带坐标）
  L1_prodigal_/L1_gms2_              → 查 02_db/<sp>/L1_coords.tsv
  L3|ext|<pid>|up<NN>aa              → 注释 CDS 向上游延伸 NN 个氨基酸
"""
import argparse
import csv
import os
import pickle
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pandas as pd

from orfscape.core.novel import iln, load_spectra, read_faa
from orfscape.core.rnaabund import (load_samples, abundance_for_span,
                                  assembled_models, tx_support)


def read_gff_cds(path):
    """locus_tag / protein_id → (contig, start, end, strand)"""
    out = {}
    for line in open(path, encoding="utf-8", errors="replace"):
        if line.startswith("#"):
            continue
        c = line.rstrip("\n").split("\t")
        if len(c) < 9 or c[2] != "CDS":
            continue
        at = {}
        for kv in c[8].split(";"):
            kv = kv.strip()
            if not kv:
                continue
            k, _, v = kv.partition(" ")
            at[k] = v.strip().strip('"')
        span = (c[0], int(c[3]), int(c[4]), c[6])
        for key in ("protein_id", "locus_tag", "ID"):
            if at.get(key):
                out[at[key]] = span
    return out


def contig_lens(fna):
    out, name, n = {}, None, 0
    for line in open(fna, encoding="utf-8", errors="replace"):
        if line.startswith(">"):
            if name:
                out[name] = n
            name, n = line[1:].split()[0], 0
        else:
            n += len(line.strip())
    if name:
        out[name] = n
    return out


class OrfMap:
    def __init__(self, root, sp="TP"):
        db = os.path.join(root, "02_db", sp)
        gen = os.path.join(root, "01_genome", sp)
        self.clen = contig_lens(os.path.join(gen, "genomic.fna"))
        self.l2len = {h.split()[0]: len(s) for h, s in
                      read_faa(os.path.join(db, "L2_6frame.faa"))}
        self.l1 = {}
        p = os.path.join(db, "L1_coords.tsv")
        if os.path.exists(p):
            for line in open(p, encoding="utf-8-sig"):
                c = line.rstrip("\n").split("\t")
                if len(c) >= 5 and c[0] != "orf":
                    try:
                        self.l1[c[0]] = (c[1], int(c[2]), int(c[3]), c[4])
                    except ValueError:
                        pass
        self.cds = read_gff_cds(os.path.join(gen, "genomic.gff"))

    def orf_span(self, orf):
        """→ (contig, strand, lo0, hi0, len_aa)  正链 0-based 闭区间"""
        s = str(orf)
        if s.startswith("L6|"):
            f = s.split("|")
            if len(f) < 4:
                return None
            a, b = f[3].split("-")
            return f[1], f[2], int(a) - 1, int(b) - 1, (int(b) - int(a) + 1) // 3
        if s.startswith("L3"):
            f = s.split("|")
            if len(f) < 4:
                return None
            pid = f[2]
            up = int("".join(ch for ch in f[3] if ch.isdigit()) or 0)
            hit = self.cds.get(pid)
            if not hit:
                return None
            ctg, cs, ce, st = hit
            if st == "+":
                lo = cs - 1 - up * 3
                hi = ce - 1
            else:
                lo = cs - 1
                hi = ce - 1 + up * 3
            return ctg, st, lo, hi, (hi - lo + 1) // 3
        if s.startswith("L1_"):
            hit = self.l1.get(s)
            return None if not hit else (hit[0], hit[3], hit[1] - 1, hit[2] - 1,
                                         (hit[2] - hit[1] + 1) // 3)
        # L2_<contig>_<frame><off>_<start>
        parts = s.split("_")
        if len(parts) < 4:
            return None
        frame = parts[-2]
        try:
            start = int(parts[-1])
        except ValueError:
            return None
        ctg = "_".join(parts[1:-2])
        L = self.l2len.get(s, 0)
        if not L:
            return None
        n3 = L * 3
        cl = self.clen.get(ctg, 0)
        if frame.startswith("+"):
            lo = start - 1
            return ctg, "+", lo, lo + n3 - 1, L
        rc_lo = start - 1
        rc_hi = rc_lo + n3 - 1
        return ctg, "-", cl - 1 - rc_hi, cl - 1 - rc_lo, L

    def pep_span(self, orf, pos1, plen):
        """肽段基因组区间（正链 0-based 闭区间）"""
        info = self.orf_span(orf)
        if not info:
            return None
        ctg, st, lo, hi, L = info
        # pFind 的 Positions 是 **0-based** 氨基酸起点（实测：ORF 首残基的肽记 0），
        # 所以偏移必须是 pos1*3；旧写法 (pos1-1)*3 会让所有肽条左移一个密码子。
        off = int(pos1) * 3
        n3 = plen * 3
        if st == "+":
            a = lo + off
            return ctg, "+", a, a + n3 - 1
        b = hi - off
        return ctg, "-", b - n3 + 1, b


def main(argv=None):
    ap = argparse.ArgumentParser(description="PT13 PDF 输入准备")
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--predict", nargs="*", default=[])
    ap.add_argument("--all-keep", action="store_true",
                    help="谱册包含全部 KEEP 行（默认只取唯一肽）")
    ap.add_argument("--rna-gtf", nargs="*", default=[],
                    help="StringTie 单样本 GTF（文件/目录/glob），用于算事件 ORF 的转录本丰度")
    args = ap.parse_args(argv)

    root = os.path.realpath(args.workdir)
    out = os.path.realpath(args.out)
    os.makedirs(out, exist_ok=True)
    om = OrfMap(root)

    # ---------- 读 P7 + 每肽最优 PSM
    v = pd.read_csv(os.path.join(root, "06b_validate_research", "novel_validated.tsv"),
                    sep="\t", dtype=str, encoding="utf-8-sig")
    v["pep"] = v["peptide"].map(iln)
    keep = v[v["final_keep"].astype(str).str.startswith("KEEP")].copy()

    best = {}
    files = []
    for sub in ("04_search_disc", "05b_research"):
        base = os.path.join(root, sub, "TP")
        if not os.path.isdir(base):
            continue
        for t in sorted(os.listdir(base)):
            f = os.path.join(base, t, "pFind-Filtered.spectra")
            if not os.path.exists(f):
                continue
            files.append(f)
            enz = "LL" if "_LL_" in t else "TL"
            for r in load_spectra(f):
                p = iln(r.get("Sequence") or "")
                if not p:
                    continue
                # 必须带 Positions（肽在 ORF 内的位点），否则画不到基因组上；
                # 以前只看 Final_Score，最优谱恰好缺 Positions 时整条肽就被丢掉
                # （TP 实测丢 4 条肽 / 2 个 ORF，事件册少 2 页）。
                if not str(r.get("Positions") or "").strip():
                    continue
                try:
                    fs = float(r.get("Final_Score") or 9e9)
                except ValueError:
                    fs = 9e9
                cur = best.get(p)
                if cur is None or fs < cur["_fs"]:
                    d = dict(r)
                    d["_fs"] = fs
                    d["_enz"] = enz
                    best[p] = d

    # ---------- Table1
    rows = []
    seen = set()
    # 事件表（P13 每肽）先读进来：Table1 的 orf 必须与事件页的 orf 一致，
    # 否则 IGV 页查不到该 ORF 的肽段（本轮踩过：肽段行画不出来）。
    p13 = os.path.join(out, "P13_3_events_per_peptide.tsv")
    if not os.path.exists(p13):
        p13 = os.path.join(os.path.dirname(out), "P13_3_events_per_peptide.tsv")
    m_ev = {}
    if os.path.exists(p13):
        _ev = pd.read_csv(p13, sep="\t", dtype=str, encoding="utf-8-sig")
        for _, r in _ev.iterrows():
            m_ev[iln(r["peptide"])] = r
    for _, r in keep.iterrows():
        p = r["pep"]
        if not args.all_keep and p in seen:
            continue
        seen.add(p)
        psm = best.get(p, {})
        e = m_ev.get(p)
        orf = str(e.get("event_orf") or "") if e is not None else ""
        if not orf:
            orf = r["protein"]
        pos = str(psm.get("Positions") or "")
        pos1 = None
        for tok in pos.replace("/", ";").split(";"):
            tok = tok.strip()
            if tok:
                try:
                    pos1 = int(tok.split(",")[0])
                except ValueError:
                    continue
                break
        plen = len(r["peptide"])
        # pos1 可以是 0（肽从 ORF 第 1 个残基起），不能用真假判断
        span = om.pep_span(orf, pos1, plen) if pos1 is not None else None
        ospan = om.orf_span(orf)
        if span and ospan:
            ctg, st, a, b = span
            rows.append(dict(
                peptide=r["peptide"], modification=(psm.get("Modification") or ""),
                orf=orf, strand=st,
                orf_start=ospan[2] + 1, orf_end=ospan[3] + 1,
                pep_nt_start=a + 1, pep_nt_end=b + 1,
                aa_offset=(pos1 + 1), mapped=1,
                tier=r.get("tier_final") or r.get("tier", ""),
                locus=r.get("locus", ""),
                n_psm=r.get("psm", ""),
                enzymes=r.get("enzymes_str", ""),
                event_type=(psm.get("_enz") and "") or "",
                event_sub="", seqid=ctg))
    t1 = pd.DataFrame(rows)
    # 事件类型：从上面已读入的 m_ev 带回来
    if m_ev:
        t1["event_type"] = t1["peptide"].map(
            lambda p: str(m_ev.get(iln(p), {}).get("event_type", "")))
        t1["event_sub"] = t1["peptide"].map(
            lambda p: str(m_ev.get(iln(p), {}).get("event_sub", "")))
    t1.to_csv(os.path.join(out, "Table1_TP_KEEP.tsv"), sep="\t", index=False,
              encoding="utf-8-sig")

    # ---------- Table3（按事件 ORF 归并）
    evrows = []
    grp = defaultdict(list)
    for r in rows:
        e = m_ev.get(iln(r["peptide"]))
        orf = e.get("event_orf") if e is not None else r["orf"]
        grp[str(orf)].append(r)
    for i, (orf, items) in enumerate(sorted(grp.items()), 1):
        info = om.orf_span(orf)
        if not info:
            continue
        ctg, st, lo, hi, L = info
        e0 = m_ev.get(iln(items[0]["peptide"]))
        evrows.append(dict(
            event_id="E%03d" % i, event_type=(e0.get("event_type") if e0 is not None else ""),
            event_subtype=(e0.get("event_sub") if e0 is not None else ""),
            orf=orf, strand=st, orf_start=lo + 1, orf_end=hi + 1,
            n_peptide=len({it["peptide"] for it in items}),
            n_psm=sum(int(it["n_psm"] or 0) for it in items),
            tier_mix="/".join(sorted({it["tier"] for it in items if it["tier"]})),
            db_mix="", seqid=ctg))
    t3 = pd.DataFrame(evrows)

    # ---------- 转录本丰度（可选）：每个事件 ORF 的位点平均 TPM
    if args.rna_gtf:
        samples = load_samples(args.rna_gtf)
        print("RNA-seq 样本 %d 个：%s" % (len(samples), ",".join(n for n, _ in samples)))
        cov, tpm, ns, txs, used = [], [], [], [], []
        for _, r in t3.iterrows():
            a = abundance_for_span(samples, r["seqid"], r["strand"],
                                   int(r["orf_start"]), int(r["orf_end"]))
            cov.append(round(a["cov"], 2) if a["cov"] == a["cov"] else "")
            tpm.append(round(a["tpm"], 1) if a["tpm"] == a["tpm"] else "")
            ns.append(a["n_sample"])
            txs.append(";".join(a["tx"]))
            used.append(a["strand_used"] or "")
        # ---- 转录本支持（证据式判据：模型是否覆盖该 ORF / 其延伸段）
        mdls = assembled_models(args.rna_gtf, log=print)

        def _tis_span(orf, strand, s0, e0):
            """与该 ORF 重叠最多的注释 CDS → 延伸段区间（无则为 None）。"""
            info = om.orf_span(orf)
            if not info:
                return None
            best = None
            for pid, (ctg, cs, ce, st) in om.cds.items():
                if st != strand or ctg != info[0]:
                    continue
                ov = min(e0, ce) - max(s0, cs) + 1
                if ov <= 0:
                    continue
                if best is None or ov > best[0]:
                    best = (ov, cs, ce)
            if not best:
                return None
            _ov, cs, ce = best
            if strand == "+":
                return (s0, cs - 1) if cs - 1 >= s0 else None
            return (ce + 1, e0) if ce + 1 <= e0 else None

        n_m, n_s, mx, ext = [], [], [], []
        for _, r in t3.iterrows():
            sup = tx_support(mdls, r["seqid"], r["strand"],
                             int(r["orf_start"]), int(r["orf_end"]))
            n_m.append(sup["n"])
            n_s.append(sup["samples"])
            mx.append(round(sup["max_tpm"], 1))
            if r["event_type"] == "annotated_correction":
                sp = _tis_span(r["orf"], r["strand"],
                               int(r["orf_start"]), int(r["orf_end"]))
                if sp:
                    e_sup = tx_support(mdls, r["seqid"], r["strand"], sp[0], sp[1])
                    ext.append("yes" if e_sup["covered"] else "no")
                else:
                    ext.append("na")
            else:
                ext.append("")
        t3["tx_models"] = n_m
        t3["tx_samples"] = n_s
        t3["tx_max_tpm"] = mx
        t3["tx_ext"] = ext
        t3["rna_cov"] = cov
        t3["rna_tpm"] = tpm
        t3["rna_n_sample"] = ns
        t3["rna_strand"] = used
        t3["rna_tx"] = txs
        t3["rna_n_total"] = len(samples)
        n_cov = sum(1 for x in cov if x != "")
        n_ok = sum(1 for x in tpm if x != "")
        print("有 RNA 覆盖的事件 ORF: %d / %d" % (n_cov, len(t3)))
    t3.to_csv(os.path.join(out, "Table3_TP_events.tsv"), sep="\t", index=False,
              encoding="utf-8-sig")

    # ---------- 谱表（合并，供引擎按 Sequence 取最优谱）
    want = {iln(r["peptide"]) for r in rows}
    out_rows = []
    for f in files:
        for r in load_spectra(f):
            if iln(r.get("Sequence") or "") in want:
                out_rows.append(r)
    if out_rows:
        cols = list(out_rows[0].keys())
        with open(os.path.join(out, "spectra_TP.tsv"), "w", newline="",
                  encoding="utf-8-sig") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t",
                               extrasaction="ignore")
            w.writeheader()
            w.writerows(out_rows)

    print("Table1 %d 行 | Table3 %d 事件 | 谱 %d 行"
          % (len(t1), len(t3), len(out_rows)))

    # ---------- 峰缓存（{(stem, scan): [(mz,int)]}，供 PDF 引擎 pickle.load 读）
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from tp_p13_tables import read_mgf_peaks, split_title
    wanted = set()
    for r in out_rows:
        stem, scan = split_title(r.get("File_Name", ""))
        wanted.add((stem, scan))
    mdir = os.path.join(root, "03_mgf")
    mgfs = sorted(os.path.join(mdir, f) for f in os.listdir(mdir)
                  if f.endswith(".mgf")) if os.path.isdir(mdir) else []
    cache = read_mgf_peaks(mgfs, wanted)
    pd.to_pickle(cache, os.path.join(out, "peaks_cache.pkl"))
    print("峰缓存 %d / %d 张谱" % (len(cache), len(wanted)))
    print("→ %s" % out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
