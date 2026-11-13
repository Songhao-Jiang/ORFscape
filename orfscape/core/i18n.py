# -*- coding: utf-8 -*-
"""ORFscape i18n：全流程与结果默认英文，可切换中文。

规则（用户 2026-11-01 指令）：
- 默认语言 = English（en）——所有日志、报告、图表、表格、Methods 均英文输出；
- 可切换中文（zh）——通过 CLI --lang zh、环境变量 ORFSCAPE_LANG=zh、或项目配置持久化。

实现方式：代码内部消息保留中文原文（作者工作语言），输出层经 L() 翻译。
  lang="en"（默认）→ 查 ZH_EN 对照表返回英文；未命中则原样输出（保底不崩）。
  lang="zh"        → 原样输出中文。
报告层（report.py）自带 EN/ZH 双标签字典，按 get_lang() 选择并切换字体。
"""
import os

_DEFAULT = "en"
_LANG = None  # 懒加载缓存

# ---------------------------------------------------------------- 消息对照表
# 键 = 代码中的中文消息原文；值 = 英文输出。
ZH_EN = {
    # 通用横幅
    "数据确定开始": "data scan started",
    "环境检测开始": "environment check started",
    "基因组获取开始": "genome fetch started",
    "电压拆分开始（D0）": "FAIMS CV split started (D0)",
    "建库开始（C 阶段）": "database build started (stage C)",
    "E1 已注释库搜库开始": "E1 annotated-DB search started",
    "QC 开始（E1 后）": "QC started (after E1)",
    "发现搜库开始": "E2 discovery search started",
    # scan
    "未发现 .raw 文件": "no .raw files found",
    "个 raw 文件": "raw file(s) discovered",
    "个文件未能推断分组": "file(s) failed pattern inference",
    "未匹配: ": "unmatched: ",
    "分组推断完成（馏分制）": "group inference done (fraction semantics)",
    "镜像酶切模式": "mirror-digest mode",
    "单酶模式": "single-enzyme mode",
    "个样品级镜像组": "sample-level mirror group(s)",
    "个预拆分电压 raw": "pre-split CV raw file(s)",
    "个多 CV 原 raw 标记跳过": "multi-CV original raw(s) marked superseded",
    "物种": "species",
    "计算 sha256（大文件较慢，可用 --no-sha256 跳过）…": "computing sha256 (slow for large files; use --no-sha256 to skip)...",
    "方法已解析": "method parsed",
    "扫描表已解析": "scan table parsed",
    "摘要": "summary",
    "馏分": "fraction(s)",
    "完成。产物目录: ": "done. output: ",
    # genome
    "获取参考装配…": "fetching reference assembly...",
    # split
    "对账基准": "reconciliation baseline",
    "待拆分 raw": "raw file(s) to split",
    "对账通过": "reconciled OK",
    "失配（见 split_report.tsv）": "mismatch(es) — see split_report.tsv",
    "拆分完成": "CV split done",
    "文件不存在: ": "file not found: ",
    "跳过含穿越成分路径: ": "skipped path with traversal component: ",
    # dbbuild
    "构建层级库…": "building layered database...",
    "声明: ": "notice: ",
    "建库完成。统计表: ": "database build done. stats: ",
    # e1/e2
    "未检测到": "not detected",
    "任务数": "task count",
    "并行模式": "parallel mode",
    "总线程": "total threads",
    "E1 完成": "E1 done",
    "E2 完成": "E2 done",
    "E2 发现搜库开始（库: ": "E2 discovery search started (dbs: ",
    # qc
    "QC（紧跟 E1）：逐任务酶切质控 + 方法/扫描警告 → 04a_qc/": "QC (right after E1): per-task cleavage QC + method/scan warnings → 04a_qc/",
    "QC（紧跟 E1）：酶切质控+门槛警告": "QC (right after E1): cleavage QC + threshold warnings",
    "QC 开始": "QC started",
    "QC 完成": "QC done",
    "条警告": "warning(s)",
    "产物 04a_qc/": "output in 04a_qc/",
    # doctor
    "检测完成": "environment check done",
    "可用": "available",
    "必需缺失": "missing required",
    "建议缺失": "missing recommended",
    "（阻断后续）": "(blocks downstream)",
    "（按降级链继续）": "(continuing via degrade chain)",
    "报告: ": "report: ",
    # novel 管线
    "消化集": "digest peptides",
    "实测靶肽": "E1 observed target peptides",
    "层命中 PSM 候选": "layer-hit PSM candidates",
    "novel 肽": "novel peptides",
    "含 decoy": "incl. decoy",
    "类内 decoy 占比": "within-class decoy fraction",
    # 检测器类目与降级说明
    "发现搜库": "discovery search",
    "谱图转换": "spectrum conversion",
    "电压拆分/转换": "CV split / conversion",
    "免费引擎(仲裁槽)": "free engine (arbitration slot)",
    "谱图预测器": "spectrum predictor",
    "功能注释": "functional annotation",
    "ab initio 预测(L1)": "ab initio prediction (L1)",
    "转录本层(L6)": "transcript layer (L6)",
    "竞争仲裁后端": "competition-arbitration backend",
    "谱图转换回退": "spectrum conversion fallback",
    "QC方法提取": "QC method extraction",
    "QC扫描表": "QC scan table",
    "仲裁引擎槽减一": "one fewer arbitration engine slot",
    "L1 层降级：无 Prodigal 则跳过 ab initio 层（六框层部分覆盖）":
        "L1 degrade: without Prodigal the ab initio layer is skipped (six-frame layer partially covers)",
    "L1 层降级：仅 Prodigal（如 Prodigal 也在缺则 L1 整层跳过）":
        "L1 degrade: Prodigal only (if Prodigal also missing, L1 is skipped entirely)",
    "L1 层降级：仅 Prodigal+Augustus": "L1 degrade: Prodigal+Augustus only",
    "L6 转录本层跳过（报告声明）": "L6 transcript layer skipped (declared in report)",
    "预测器回退 pDeep2（如配）或无预测模式": "predictor falls back to pDeep2 (if configured) or no-prediction mode",
    "功能注释降级（跳过 BLASTp 层）": "functional annotation degraded (BLASTp layer skipped)",
    # QC 摘要
    "QC 摘要（E1 后）": "QC summary (after E1)",
    "任务数（raw×CV）": "tasks (raw x CV)",
    "警告数": "warning count",
    "警告": "warnings",
    "无": "none",
    "酶切特异性(两端全切近似)": "cleavage specificity (approx. fully-cleaved)",
}


def get_lang():
    """当前输出语言（默认 en）。优先级：set_lang 缓存 > ORFSCAPE_LANG 环境变量 > en。"""
    global _LANG
    if _LANG is None:
        _LANG = os.environ.get("ORFSCAPE_LANG", _DEFAULT).strip().lower() or _DEFAULT
    return _LANG


def set_lang(lang):
    """切换语言（'en'/'zh'），同步写 ORFSCAPE_LANG 供子进程继承。"""
    global _LANG
    lang = (lang or _DEFAULT).strip().lower()
    if lang not in ("en", "zh"):
        lang = _DEFAULT
    _LANG = lang
    os.environ["ORFSCAPE_LANG"] = lang
    return lang


def _zh_keys_longest_first():
    return sorted(ZH_EN.keys(), key=len, reverse=True)


_KEYS_CACHE = None


def L(text):
    """输出层翻译（子串替换式，兼容 f-string 内嵌中文）。

    lang=en（默认）：把文本中出现的每个中文片段按键长优先替换为英文；
    未命中片段保持原样（保底不崩）。lang=zh：原样返回。
    """
    if get_lang() != "en":
        return text
    global _KEYS_CACHE
    if _KEYS_CACHE is None:
        _KEYS_CACHE = _zh_keys_longest_first()
    out = text
    for k in _KEYS_CACHE:
        if k in out:
            out = out.replace(k, ZH_EN[k])
    return out


def lang_name():
    return {"en": "English", "zh": "中文"}[get_lang()]
