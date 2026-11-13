# -*- coding: utf-8 -*-
"""真核示例物种的"文档参数"预设。

数字与流程全部来自毕业设计材料，不在代码里臆造：

- **CREI 莱茵衣藻**：`材料/成对标记/项目全流程文档_20261111.md` §2.2
  CC-503 v5.5（GCF_000002595.2）、17 条染色体、**19,527 注释蛋白**、
  部分标记 0|100 数据集 PRIDE PXD003236、权重/阈值**零物种特调**迁移验证。
- **HSAP 人**：`材料/人睾丸蛋白质基因组/人睾丸暗蛋白质组_全流程实验设计_v2.0_20261031.md` §一
  GCF_000001405.40（RS_2025_08）、ref1（UP000005640 一基因一蛋白）20,652、
  db0（RefSeq 136,807 → I/L 去重 90,373 + CON 571）90,944、
  6,384,275 转录本 → 2,272,760 高置信 ORF、人睾丸 36 Trypsin + 36 LysargiNase MGF。

预设有两个用途：
1. `orfscape eukbuild --preset CREI` 自动定位 `01_genome/<code>/` 三件套；
2. 建完自动核对**文档里的蛋白数**，偏差 >1% 就在 manifest 与日志里标出来。
"""
PRESETS = {
    "CREI": {
        "code": "CREI",
        "name": "Chlamydomonas reinhardtii",
        "assembly": "GCF_000002595.2_Chlamydomonas_reinhardtii_v5.5",
        "expected_proteins": 19527,
        "chromosomes": 17,
        "enzymes": ["TL", "LL"],
        "dataset": "PRIDE PXD003236（0|100 部分标记，cw15）",
        "source": "成对标记/项目全流程文档_20261111.md §2.2",
        "notes": "权重与阈值零物种特调迁移；576 候选谱 → 179 同位素支持 → 39 事件",
    },
    "HSAP": {
        "code": "HSAP",
        "name": "Homo sapiens",
        "assembly": "GCF_000001405.40_GRCh38.p14",
        "expected_proteins": None,
        "ref1_proteins": 20652,
        "db0_proteins": 90944,
        "transcripts": 6384275,
        "high_conf_orf": 2272760,
        "enzymes": ["TL", "LL"],
        "dataset": "人睾丸 36 Trypsin（48 GB）+ 36 LysargiNase（37 GB）MGF",
        "source": "人睾丸蛋白质基因组/人睾丸暗蛋白质组_全流程实验设计_v2.0_20261031.md §一",
        "notes": "ref1=UP000005640 一基因一蛋白 20,652；db0=RefSeq 136,807 → I/L 去重 90,373 + CON 571",
    },
}

ALIASES = {"CHLAMY": "CREI", "CHLAMYDOMONAS": "CREI", "HUMAN": "HSAP",
           "HOMO": "HSAP", "衣藻": "CREI", "人": "HSAP"}


def get(name):
    """按代码/别名取预设；未知返回 None。"""
    if not name:
        return None
    key = str(name).strip().upper()
    key = ALIASES.get(key, key)
    return PRESETS.get(key)


def paths(workdir, code):
    """`01_genome/<code>/` 下的三件套路径（可能不存在，由调用方校验）。"""
    import os
    base = os.path.join(workdir, "01_genome", code)
    return {"fna": os.path.join(base, "genomic.fna"),
            "gff": os.path.join(base, "genomic.gff"),
            "ref": os.path.join(base, "protein.faa")}


def check(observed, code=None, expected=None):
    """把实测蛋白数与文档值比对；返回 {expected, observed, delta_pct, ok}。"""
    preset = get(code) or {}
    want = expected if expected is not None else preset.get("expected_proteins")
    if not want:
        return {"expected": None, "observed": observed, "delta_pct": None,
                "ok": True}
    delta = (observed - want) / float(want) * 100.0
    return {"expected": want, "observed": observed,
            "delta_pct": round(delta, 3), "ok": abs(delta) <= 1.0}
