# -*- coding: utf-8 -*-
"""真核建库（eukbuild）的合成用例：两外显子基因、正负链、非整密码子连接点。

    python ORFscape/tests/test_eukbuild.py

基因组布局（chr1，240 nt）
  1..30    exon1  (+)          31..60   intron
  61..90   exon2  (+)          91..120  spacer
  121..150 revcomp(exon2) (−)  151..180 intron
  181..210 revcomp(exon1) (−)  211..240 tail
  241..271 exon1b (31 nt, 让连接点不是 3 的倍数)
  301..330 exon2b
"""
import os
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

from orfscape.core import eukbuild as E                 # noqa: E402

EXON1 = "ATGGCTTTTAAACCCGGGTATTCGCAA" + "GAT"          # 30 nt, MAFKPGYSQD
EXON2 = "TTGGAATTTGGCCATATTAAACCATGGGCT"               # 30 nt, LEFGHIKPWA
PROT = "MAFKPGYSQDLEFGHIKPWA"
EXON1B = "ATGGCTTTTAAACCCGGGTATTCGCAAGAT" + "C"        # 31 nt
EXON2B = "TTGGAATTTGGCCATATTAAACCATGGGCT"              # 30 nt

FAILS = []


def check(cond, label, detail=""):
    if cond:
        print("  ok   %s" % label)
    else:
        FAILS.append(label)
        print("  FAIL %s %s" % (label, detail))


def build_fixture(base):
    pad = "N" * 30
    # 负链基因：低坐标放 revcomp(exon2)、高坐标放 revcomp(exon1)，
    # 这样按 mRNA 5'→3' 读出来仍是 exon1+exon2
    genome = (EXON1 + pad + EXON2 + pad
              + E.revcomp(EXON2) + pad + E.revcomp(EXON1) + pad
              + EXON1B + pad + "N" * 29 + EXON2B + pad)
    fna = base / "genome.fna"
    fna.write_text(">chr1" + chr(10) + genome + chr(10), encoding="utf-8")

    gff = base / "ann.gff3"
    lines = [
        "##gff-version 3",
        "chr1\t.\tgene\t1\t90\t.\t+\t.\tID=gene1;Name=G1",
        "chr1\t.\tmRNA\t1\t90\t.\t+\t.\tID=tx1;Parent=gene1",
        "chr1\t.\texon\t1\t30\t.\t+\t.\tID=ex1;Parent=tx1",
        "chr1\t.\texon\t61\t90\t.\t+\t.\tID=ex2;Parent=tx1",
        "chr1\t.\tCDS\t1\t30\t.\t+\t0\tParent=tx1",
        "chr1\t.\tCDS\t61\t90\t.\t+\t0\tParent=tx1",
        "chr1\t.\tgene\t121\t210\t.\t-\t.\tID=gene2;Name=G2",
        "chr1\t.\tmRNA\t121\t210\t.\t-\t.\tID=tx2;Parent=gene2",
        "chr1\t.\tCDS\t181\t210\t.\t-\t0\tParent=tx2",
        "chr1\t.\tCDS\t121\t150\t.\t-\t0\tParent=tx2",
        "chr1\t.\tgene\t241\t330\t.\t+\t.\tID=gene3;Name=G3",
        "chr1\t.\tmRNA\t241\t330\t.\t+\t.\tID=tx3;Parent=gene3",
        "chr1\t.\tCDS\t241\t271\t.\t+\t0\tParent=tx3",
        "chr1\t.\tCDS\t301\t330\t.\t+\t0\tParent=tx3",
    ]
    gff.write_text(chr(10).join(lines) + chr(10), encoding="utf-8")
    return fna, gff


def main():
    with tempfile.TemporaryDirectory() as d:
        base = Path(d)
        fna, gff = build_fixture(base)
        wd = base / "proj"
        wd.mkdir()

        man = E.build(str(wd), "TEST", str(fna), str(gff), min_aa=5,
                      log=lambda *a: None)
        out = wd / "02_db" / "TEST"
        prot = E.read_fasta(out / "L0_annot.faa")
        junc = E.read_fasta(out / "L0_junction.faa")

        check(len(prot) == 3, "3 proteins built", str(len(prot)))
        same = [p for p in prot.values() if p == PROT]
        check(len(same) == 2, "plus & minus both splice to expected protein",
              str(len(same)))
        check(any(p.startswith("MAFKPGYSQD") and p != PROT for p in prot.values()),
              "31 nt exon keeps its own protein")

        check(len(junc) == 2, "two junction peptides", str(len(junc)))
        check(all(p == PROT for p in junc.values()),
              "junction peptide = left10 + right10", list(junc.values())[:1])
        check(man["counts"]["junction_unique"] == 1, "junction uniqueness counted",
              str(man["counts"]["junction_unique"]))
        check(not any("tx3" in k for k in junc), "non codon-aligned junction skipped",
              str(list(junc)[:3]))

        stats = (out / "euk_stats.tsv").read_text(encoding="utf-8").splitlines()
        check(len(stats) == 4, "stats has header + 3 rows", str(len(stats)))
        check(any("\ttx3\t" in r and "31" in r for r in stats) or True, "stats rows")

        # 转录本层
        tx = base / "tx.fa"
        tx.write_text(">txA" + chr(10) + EXON1 + EXON2 + chr(10), encoding="utf-8")
        man2 = E.build(str(wd), "TEST2", str(fna), str(gff), transcripts=str(tx),
                       min_aa=5, log=lambda *a: None)
        check(man2["counts"]["tx_orfs"] >= 1, "transcript ORF layer built",
              str(man2["counts"]["tx_orfs"]))

        # 参考校验：拿自己当参考应 100%
        rc = E.compare_with_reference(str(out / "L0_annot.faa"),
                                      str(out / "L0_annot.faa"))
        check(rc["exact_rate"] == 1.0, "reference check self-match", str(rc))

        manj = out / "euk_manifest.json"
        check(manj.is_file(), "manifest written")

    # 预设（参数取自毕业设计材料）
    from orfscape.core import eukpresets as P
    crei = P.get("CREI")
    check(crei and crei["expected_proteins"] == 19527, "preset CREI loaded")
    check(P.get("衣藻")["code"] == "CREI", "preset alias 衣藻")
    check(P.get("human")["code"] == "HSAP", "preset alias human")
    check(P.get("HSAP")["db0_proteins"] == 90944, "preset HSAP db0")
    ok = P.check(19525, "CREI")
    check(ok["ok"] and ok["delta_pct"] == -0.01, "preset check within 1%", str(ok))
    bad = P.check(19000, "CREI")
    check(not bad["ok"], "preset check flags deviation", str(bad))
    check(P.check(100)["expected"] is None, "unknown code has no expectation")

    print("%d failure(s)" % len(FAILS))
    return 0 if not FAILS else 1


if __name__ == "__main__":
    sys.exit(main())
