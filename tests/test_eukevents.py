# -*- coding: utf-8 -*-
"""真核事件分类（eukevents）的合成用例：七类判定 + 跨外显子 + 肽级归一。

    python ORFscape/tests/test_eukevents.py

合成模型（NC_1，正链转录本 X1，200 nt 拼接长度）
  exon1 101..200 / exon2 301..400      CDS 131..200 + 301..380（两端各留 UTR）
  连接点 nt=100（exon1 末尾）          → LJR 肽跨 165..200 + 301..332
"""
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

from orfscape.core import eukevents as E                 # noqa: E402

FAILS = []


def check(cond, label, detail=""):
    if cond:
        print("  ok   %s" % label)
    else:
        FAILS.append(label)
        print("  FAIL %s %s" % (label, detail))


def tx(tid_seq="NC_1", strand="+", exons=None, cds=None):
    exons = exons or [(101, 200), (301, 400)]
    cds = cds if cds is not None else [(131, 200, 0), (301, 380, 0)]
    return {"seqid": tid_seq, "strand": strand, "gene": "G1", "biotype": "mRNA",
            "exons": list(exons), "cds": list(cds),
            "nt_len": sum(e - s + 1 for s, e in exons)}


def gene(name, seqid="NC_1", strand="+", lo=101, hi=400, cds=None):
    cds = cds if cds is not None else [(131, 200, 0), (301, 380, 0)]
    exons = [(lo, hi)] if lo <= 200 else []
    return {"seqid": seqid, "strand": strand, "biotype": "mRNA", "gene": name,
            "exons": exons or [(101, 200), (301, 400)], "cds": cds,
            "lo": lo, "hi": hi}


def genes(**items):
    from collections import defaultdict
    out = defaultdict(list)
    for g in items.values():
        out[g["seqid"]].append(g)
    return out


def main():
    X1 = tx()
    txs = {"X1": X1}

    # ---- 1. 头部解析
    kind, tid, extra = E.parse_euk_name(
        "L2|X1|f2|nt=33|len=80")
    check(kind == "L2" and tid == "X1" and extra.get("f") == "f2", "parse L2 header",
          str((kind, tid, extra)))
    check(E.parse_euk_name("LJR|X1|nt=120|len=24")[0] == "LJR", "parse LJR header")
    check(E.parse_euk_name("L0|X1|gene=G1|chr=NC_1")[2].get("gene") == "G1",
          "parse L0 header")
    check(E.parse_euk_name("TP_001") is None, "reject non-euk header")

    # ---- 2. 坐标映射（正链）
    b = E.tx_blocks(X1, 0, 9)
    check(b == [("NC_1", "+", 100, 108)], "plus nt0-9 -> 101..109 (0-based)", str(b))
    b = E.tx_blocks(X1, 195, 6)                 # exon2 内偏移 95..100（末尾截断到 200）
    check(b == [("NC_1", "+", 395, 399)], "plus nt195-200 -> 396..400", str(b))

    # ---- 3. 跨外显子（LJR 肽）
    b = E.pep_blocks("LJR|X1|nt=100|len=24", 0, 24, txs)
    check(len(b) == 2, "junction peptide spans two blocks", str(b))
    check(b[0][2:] == (164, 199) and b[1][2:] == (300, 335),
          "junction block coordinates", str(b))

    # ---- 4. 分类：注释内同框
    cls, g = E.classify_peptide("L2|X1|f1|nt=100|len=25", 10, 10, txs, genes(G1=gene("G1")))
    check(cls == E.INTERNAL and g == "G1", "in-frame inside CDS -> internal", cls)

    # ---- 5. N 端延伸（肽段在 CDS 上游 UTR 内）
    cls, _g = E.classify_peptide("L2|X1|f1|nt=0|len=80", 0, 3, txs, genes(G1=gene("G1")))
    check(cls == E.NTERM, "peptide upstream of CDS -> nterm_ext", cls)

    # ---- 6. C 端延伸（肽段越过 CDS 末端）
    cls, _g = E.classify_peptide("L2|X1|f1|nt=180|len=10", 2, 3, txs, genes(G1=gene("G1")))
    check(cls == E.CTERM, "peptide downstream of CDS -> cterm_ext", cls)

    # ---- 7. ORF 延伸但肽段在 CDS 内 → alt_frame
    cls, _g = E.classify_peptide("L2|X1|f1|nt=0|len=80", 45, 3, txs, genes(G1=gene("G1")))
    check(cls == E.ALT_FRAME, "ORF extends, peptide inside -> alt_frame", cls)

    # ---- 8. 基因间区（该 contig 无基因）
    X2 = dict(X1, seqid="NC_2", gene="")
    cls, _g = E.classify_peptide("L0|X2|gene=|chr=NC_2", 10, 3, {"X2": X2},
                                 genes(G1=gene("G1")))
    check(cls == E.INTERGENIC, "no gene -> intergenic", cls)

    # ---- 9. 反义链
    X3 = dict(X1, seqid="NC_3")
    G3 = gene("G3", seqid="NC_3", strand="-")
    cls, _g = E.classify_peptide("L0|X3|gene=G3", 10, 3, {"X3": X3}, genes(G3=G3))
    check(cls == E.ANTISENSE, "opposite strand gene -> antisense", cls)

    # ---- 10. ncRNA（无 CDS 的基因）
    X4 = dict(X1, seqid="NC_4")
    G4 = gene("G4", seqid="NC_4", cds=[])
    cls, _g = E.classify_peptide("L2|X4|f1|nt=0|len=80", 45, 3, {"X4": X4}, genes(G4=G4))
    check(cls == E.NCRNA, "gene without CDS -> ncRNA_translation", cls)

    # ---- 11. Ig/TCR 基因座排除
    X5 = dict(X1, seqid="NC_5")
    G5 = gene("IGHG1", seqid="NC_5")
    cls, _g = E.classify_peptide("L2|X5|f1|nt=100|len=25", 10, 3, {"X5": X5}, genes(G5=G5))
    check(cls == E.INTERGENIC, "Ig/TCR gene excluded", cls)

    # ---- 12. 肽级归一：Tier 优先，其次有基因归属
    hits = {"PEP": [
        {"class": E.ALT_FRAME, "gene": "", "tier": "T3", "orf": "a", "pos": 1},
        {"class": E.NTERM, "gene": "G1", "tier": "T1", "orf": "b", "pos": 2},
    ]}
    main = E.normalize_peptide_level(hits)
    check(main["PEP"]["class"] == E.NTERM, "peptide-level: best tier wins",
          str(main["PEP"]))
    check(main["PEP"]["n_events"] == 2, "peptide-level: event count kept")
    hits2 = {"P2": [
        {"class": E.INTERGENIC, "gene": "", "tier": "T1", "orf": "a", "pos": 1},
        {"class": E.INTERNAL, "gene": "G1", "tier": "T1", "orf": "b", "pos": 2},
    ]}
    check(E.normalize_peptide_level(hits2)["P2"]["class"] == E.INTERNAL,
          "peptide-level: gene-bearing wins on tie")

    # ---- 13. 汇总只报出现过的类
    s = E.summarize(main)
    check(s.get(E.NTERM) == 1 and E.INTERGENIC not in s, "summarize counts only present")

    # ---- 14. GFF → 模型（1 转录本 / 2 外显子 / 2 CDS）
    with tempfile.TemporaryDirectory() as tmp:
        gff = Path(tmp) / "mini.gff3"
        rows = [
            "##gff-version 3",
            "NC_1\ttest\tgene\t101\t400\t.\t+\t.\tID=gene1;Name=G1;biotype=mRNA",
            "NC_1\ttest\tmRNA\t101\t400\t.\t+\t.\tID=X1;Parent=gene1",
            "NC_1\ttest\texon\t101\t200\t.\t+\t.\tID=e1;Parent=X1",
            "NC_1\ttest\texon\t301\t400\t.\t+\t.\tID=e2;Parent=X1",
            "NC_1\ttest\tCDS\t131\t200\t.\t+\t0\tID=c1;Parent=X1",
            "NC_1\ttest\tCDS\t301\t380\t.\t+\t0\tID=c2;Parent=X1",
        ]
        gff.write_text("\n".join(rows) + "\n", encoding="utf-8")
        txs2, genes2 = E.load_tx_model(str(gff))
        check(len(txs2) == 1 and txs2["X1"]["nt_len"] == 200, "GFF: one transcript, 200 nt",
              str({k: v["nt_len"] for k, v in txs2.items()}))
        check(len(txs2["X1"]["cds"]) == 2, "GFF: two CDS rows")
        g2 = genes2.get("gene1")
        check(g2 and g2["lo"] == 101 and g2["hi"] == 400, "GFF: gene span", str(g2 and (g2["lo"], g2["hi"])))

    # ---- 15. build_novgene_euk 端到端（临时项目：GFF + 层蛋白 + 谱表 + KEEP 表）
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        (root / "01_genome" / "TEST").mkdir(parents=True)
        (root / "02_db" / "TEST").mkdir(parents=True)
        (root / "04_search_disc" / "TEST" / "t1").mkdir(parents=True)
        (root / "06_validate").mkdir(parents=True)
        (root / "01_genome" / "TEST" / "genomic.gff").write_text(
            "\n".join([
                "##gff-version 3",
                "NC_1\ttest\tgene\t101\t400\t.\t+\t.\tID=gene1;Name=G1",
                "NC_1\ttest\tmRNA\t101\t400\t.\t+\t.\tID=X1;Parent=gene1",
                "NC_1\ttest\texon\t101\t200\t.\t+\t.\tID=e1;Parent=X1",
                "NC_1\ttest\texon\t301\t400\t.\t+\t.\tID=e2;Parent=X1",
                "NC_1\ttest\tCDS\t131\t200\t.\t+\t0\tID=c1;Parent=X1",
                "NC_1\ttest\tCDS\t301\t380\t.\t+\t0\tID=c2;Parent=X1",
            ]) + "\n", encoding="utf-8")
        (root / "02_db" / "TEST" / "euk_manifest.json").write_text(
            '{"species": "TEST"}\n', encoding="utf-8")
        spectra = [
            "Sequence\tProteins\tPositions\tTarget/Decoy",
            "PEPTADEAKAA\tL2|X1|f1|nt=100|len=25\t10\tTarget",
            "PEPTADB\tL2|X1|f1|nt=0|len=80\t0\tTarget",
            "PEPTADECCCC\tL0|X1|gene=G1|chr=NC_1|strand=+|exons=2|len=90\t10\tTarget",
        ]
        (root / "04_search_disc" / "TEST" / "t1" / "pFind-Filtered.spectra").write_text(
            "\n".join(spectra) + "\n", encoding="utf-8")
        (root / "06_validate" / "novel_validated.tsv").write_text(
            "\n".join([
                "species\tpeptide\tfinal_keep\ttier",
                "TEST\tPEPTADEAKAA\tKEEP\tT1",
                "TEST\tPEPTADB\tKEEP\tT2",
                "TEST\tPEPTADECCCC\tDROP\tT3",
            ]) + "\n", encoding="utf-8")
        main2 = E.build_novgene_euk(str(root), str(root / "07_novgene"), log=lambda *_: None)
        check(len(main2) == 2, "P8 euk: only KEEP rows classified", str(len(main2)))
        check(main2[("TEST", "PEPTADEAKAA")]["class"] == E.INTERNAL,
              "P8 euk: internal class", str(main2.get(("TEST", "PEPTADEAKAA"))))
        check(main2[("TEST", "PEPTADB")]["class"] == E.NTERM,
              "P8 euk: nterm class", str(main2.get(("TEST", "PEPTADB"))))
        check((root / "07_novgene" / "euk_events.tsv").is_file(), "P8 euk: events table written")
        check((root / "07_novgene" / "euk_events_summary.tsv").is_file(),
              "P8 euk: summary written")

    print("%d failure(s)" % len(FAILS))
    return 0 if not FAILS else 1


if __name__ == "__main__":
    sys.exit(main())
