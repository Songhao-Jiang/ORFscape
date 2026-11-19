# -*- coding: utf-8 -*-
"""六框适用性判定（sframe）的合成用例。

    python ORFscape/tests/test_sframe.py

构造两个合成基因组：
  compact   CDS 覆盖 60%（R≈10）→ 应判 compact
  sparse    CDS 覆盖 6%（R≈100）→ 应判 extreme
并验证六框翻译的 ORF 计数与基因间拆分。
"""
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

from orfscape.core import sframe as S                    # noqa: E402

FAILS = []


def check(cond, label, detail=""):
    if cond:
        print("  ok   %s" % label)
    else:
        FAILS.append(label)
        print("  FAIL %s %s" % (label, detail))


def make_project(tmp, name, cds_frac, cds_len=600, genome_len=1000):
    """一个 contig + 一个基因（CDS 占 cds_frac）。"""
    root = Path(tmp)
    gdir = root / "01_genome" / name
    gdir.mkdir(parents=True)
    # 无终止子的可翻译序列（全 GCT/GCG 之类），保证六框里能出长 ORF
    seq = "GCT" * (genome_len // 3)
    (gdir / "genomic.fna").write_text(">chr1\n%s\n" % seq, encoding="utf-8")
    gstart, gend = 1, int(genome_len * cds_frac)
    lines = [
        "chr1\tsyn\tgene\t%d\t%d\t.\t+\t.\tID=gene-1;gene_biotype=protein_coding"
        % (gstart, gend),
        "chr1\tsyn\tCDS\t%d\t%d\t.\t+\t0\tID=cds-1;Parent=gene-1"
        % (gstart, min(gend, cds_len)),
    ]
    (gdir / "genomic.gff").write_text("\n".join(lines) + "\n", encoding="utf-8")
    # 注释蛋白：长度 = CDS/3，全部为 A
    n_aa = max(1, (min(gend, cds_len) - gstart + 1) // 3)
    (gdir / "protein.faa").write_text(">p1\n%s\n" % ("A" * n_aa), encoding="utf-8")
    return root, name


def main():
    with tempfile.TemporaryDirectory() as tmp:
        # 1) 紧凑基因组：CDS 60% → R = 2G/P ≈ 6/0.6 = 10 → compact
        root, name = make_project(tmp, "cmp", cds_frac=0.60)
        m = S.genome_metrics(root, name)
        v = S.verdict(m)
        check(m["genome_bp"] == 999, "基因组长度按 fasta 统计", m["genome_bp"])
        check(abs(v["coding_density"] - 0.6) < 0.02, "编码密度 ≈ 60%",
              v["coding_density"])
        check(v["level"] == "compact", "紧凑基因组判 compact", v)

        # 2) 稀疏基因组：CDS 6% → R ≈ 100 → extreme
        root2, name2 = make_project(Path(tmp) / "b", "spa", cds_frac=0.06,
                                    cds_len=60, genome_len=1000)
        m2 = S.genome_metrics(root2, name2)
        v2 = S.verdict(m2)
        check(v2["level"] == "extreme", "稀疏基因组判 extreme", v2)

        # 3) 六框实测：全 GCT 的序列没有终止子 → 每条链 3 帧各 1 条长 ORF
        fa = Path(root) / "01_genome" / name / "genomic.fna"
        m6 = S.sixframe_measure(str(fa), min_aa=30)
        check(m6["sixframe_orfs"] == 6, "六框各出一条 ORF", m6["sixframe_orfs"])
        check(m6["orfs_with_x"] == 0, "无含 X 的 ORF", m6)

        # 4) 基因间拆分：基因占 1..600，正链帧 0 的 ORF 从 0 起 → 部分重叠
        gff = Path(root) / "01_genome" / name / "genomic.gff"
        m6b = S.sixframe_measure(str(fa), min_aa=30, gene_gff=str(gff))
        check("sixframe_intergenic_orfs" in m6b, "给出基因间 ORF 数", m6b)
        check(0 <= m6b["sixframe_intergenic_orfs"] <= m6b["sixframe_orfs"],
              "基因间 ORF 不超过总数", m6b)

        # 5) 算力闸：把上限压到 1 条，必须判 extreme
        old = S.ABS_ORF_MAX
        try:
            S.ABS_ORF_MAX = 1
            v3 = S.verdict(m6)
            check(v3["level"] == "extreme" and v3["blocked_by"],
                  "超过 ORF 上限时强制不做", v3)
        finally:
            S.ABS_ORF_MAX = old

        # 6) L1_SF 名字解析 + 基因组坐标映射（正负链）
        from orfscape.core import eukevents as E
        sf = E._l1sf("L1_SF|chr1|+f1|nt=100|len=50")
        check(sf == ("chr1", "+", 100, 50), "L1_SF 名字解析", sf)
        check(E._l1sf("L1_SF|chr1|x9|nt=1|len=2") is None, "非法名字返回 None")

        plus = E.pep_blocks("L1_SF|chr1|+f1|nt=100|len=50", 2, 10, {})
        check(plus == [("chr1", "+", 106, 135)], "正链肽段区间", plus)
        minus = E.pep_blocks("L1_SF|chr1|-f2|nt=1000|len=50", 2, 10, {})
        check(minus == [("chr1", "-", 965, 994)], "负链肽段区间", minus)
        ob = E.orf_blocks("L1_SF|chr1|-f2|nt=1000|len=50", {})
        check(ob == [("chr1", "-", 851, 1000)], "负链 ORF 跨度", ob)

    print("")
    if FAILS:
        print("%d failure(s): %s" % (len(FAILS), ", ".join(FAILS)))
        return 1
    print("0 failure(s)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
