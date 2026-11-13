# -*- coding: utf-8 -*-
"""ORFscape 全局配置：工具路径、物种别名表、酶表、raw 文件名分组模式库。"""
import re
from pathlib import Path

# ---------------------------------------------------------------- 工具
TOOLSDIR = Path(__file__).resolve().parent / "tools"
RAWMETHOD_EXE = TOOLSDIR / "rawmethod" / "rawmethod.exe"   # raw 内嵌方法提取（OpenMcdf，无需 Thermo 运行库）
RAWDUMP_EXE = TOOLSDIR / "rawdump" / "rawdump.exe"         # 逐扫描表导出（Thermo RawFileReader）

# ---------------------------------------------------------------- 物种别名表
# 键为文件名/用户输入中常见写法（大小写不敏感，查表前统一 lower）。
SPECIES_ALIAS = {
    "mab": "Mycobacterium abscessus",
    "m. abscessus": "Mycobacterium abscessus",
    "abscessus": "Mycobacterium abscessus",
    "脓肿分枝杆菌": "Mycobacterium abscessus",
    "mch": "Mycobacterium chelonae",
    "m. chelonae": "Mycobacterium chelonae",
    "chelonae": "Mycobacterium chelonae",
    "龟分枝杆菌": "Mycobacterium chelonae",
    # —— 后续物种在此扩展 ——
    "eco": "Escherichia coli",
    "e. coli": "Escherichia coli",
    "ecoli": "Escherichia coli",
    "mtb": "Mycobacterium tuberculosis",
    "smeg": "Mycobacterium smegmatis",
}

# ---------------------------------------------------------------- NCBI 检索别名链
# 2018 年分枝杆菌大量改名（Mycobacterium → Mycobacteroides 等），旧名在 esearch
# 可能不再命中；每个规范名给一组按序尝试的检索词，首个返回结果者胜出。
TAX_SEARCH_ALIASES = {
    "Mycobacterium abscessus": ["Mycobacterium abscessus", "Mycobacteroides abscessus"],
    "Mycobacterium chelonae": ["Mycobacterium chelonae", "Mycobacteroides chelonae"],
}

# ---------------------------------------------------------------- 酶表
# cleavage: KR_C = C 端 K/R 后切（Trypsin 类）；KR_N = N 端 K/R 前切（LysargiNase 类）
# pfind_enzyme: pFind 酶名（LyasrgiNase 为 pFind 历史遗留拼写，调度时自动处理）
ENZYMES = {
    "TL": {"name": "Trypsin+LysC", "cleavage": "KR_C", "pfind_enzyme": "Trypsin", "mirror_of": "LL"},
    "LL": {"name": "LysargiNase+LysN", "cleavage": "KR_N", "pfind_enzyme": "LysargiNase", "mirror_of": "TL"},
    "Trypsin": {"name": "Trypsin", "cleavage": "KR_C", "pfind_enzyme": "Trypsin"},
    "LysC": {"name": "LysC", "cleavage": "KR_C", "pfind_enzyme": "LysC"},
    "LysN": {"name": "LysN", "cleavage": "KR_N", "pfind_enzyme": "LysN"},
    "LysargiNase": {"name": "LysargiNase", "cleavage": "KR_N", "pfind_enzyme": "LysargiNase"},
    "GluC": {"name": "GluC", "cleavage": "DE_C", "pfind_enzyme": "GluC"},
    "AspN": {"name": "AspN", "cleavage": "D_N", "pfind_enzyme": "AspN"},
    "ArgC": {"name": "ArgC", "cleavage": "R_C", "pfind_enzyme": "ArgC"},
    "Chymotrypsin": {"name": "Chymotrypsin", "cleavage": "FWYL_C", "pfind_enzyme": "Chymotrypsin"},
}
# 镜像酶对（同一样品跨这两种酶 → 镜像组）
MIRROR_PAIRS = {frozenset({"TL", "LL"})}

# ---------------------------------------------------------------- raw 文件名分组模式库
# 按序尝试；命中即停。map_* 为捕获组取值映射；缺省 condition = sample_type。
RAW_PATTERNS = [
    {
        "name": "ntm_full",
        "rx": re.compile(r"_(?P<sp>MAB|Mab|MCH|Mch)_(?P<st>T|C)_(?P<enz>TL|LL)_M(?P<rep>\d+)\.raw$", re.IGNORECASE),
        "map_species": {"mab": "Mab", "Mab": "Mab", "MAB": "Mab", "MCH": "Mch", "Mch": "Mch"},
        "map_sample": {"T": "TCL", "C": "CFP"},
    },
    {
        "name": "ntm_mixture",
        # Mab 全菌裂解线：MAB_T_Mixture_1 / MAB_L_Mixture_1（酶为单字母 T/L）
        "rx": re.compile(r"_(?P<sp>MAB|Mab|MCH|Mch)_(?P<enzs>T|L)_Mixture_(?P<rep>\d+)\.raw$", re.IGNORECASE),
        "map_species": {"MAB": "Mab", "Mab": "Mab", "MCH": "Mch", "Mch": "Mch"},
        "map_enz": {"T": "TL", "L": "LL"},
        "sample_type": "TCL",
    },
    {
        "name": "tp_mirror",
        # 梅毒螺旋体：..._withtrap_T1_45.raw / _L1_65.raw
        # T=Trypsin, L=LysargiNase；末尾 45/65 为 FAIMS 电压；物种固定 TP（文件名不含物种）
        "rx": re.compile(r"_(?P<enzs>T|L)(?P<rep>\d+)_(?:45|65)\.raw$", re.IGNORECASE),
        "species": "TP",
        "map_enz": {"T": "TL", "L": "LL"},
        "sample_type": "TP",
    },
]

# ---------------------------------------------------------------- 目录约定
PROJECT_SUBDIRS = [
    "00_input", "01_genome", "02_db", "02_tx", "03_mgf",
    "04_search_annot", "04a_qc", "04_search_disc",
    "05_novel", "06_arbitration", "07_annotation", "08_report",
]
