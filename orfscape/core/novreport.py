# -*- coding: utf-8 -*-
"""新基因事件汇总报告（Word）：事件分类 + 结构/功能预测 + 同源比对。

数据来源（08_report 下）：
  NOV_screening_table.tsv、NG_table.tsv、AN_table.tsv
  NG1/NG2、AN1–AN4 图
输出：NovelGene_Report.docx（EN）/ NovelGene_Report_zh.docx（ZH）
"""
import datetime
import os
from collections import Counter

from .report import load_tsv


def _num(v, default=0.0):
    try:
        return float(str(v).replace(",", ""))
    except (TypeError, ValueError):
        return default


def _median(vals):
    v = sorted(x for x in vals if x is not None)
    if not v:
        return 0.0
    n = len(v)
    return v[n // 2] if n % 2 else 0.5 * (v[n // 2 - 1] + v[n // 2])


def build_novelreport(root, out_dir, log=print):
    """生成 Word 报告；语言跟随 get_lang()。"""
    from docx import Document
    from docx.enum.text import WD_ALIGN_PARAGRAPH
    from docx.oxml.ns import qn
    from docx.shared import Pt, Inches

    from .i18n import get_lang
    zh = get_lang() == "zh"

    ng = load_tsv(os.path.join(out_dir, "NG_table.tsv"))
    an = load_tsv(os.path.join(out_dir, "AN_table.tsv"))
    scr = load_tsv(os.path.join(out_dir, "NOV_screening_table.tsv"))

    # 列名归一化：中英表头 → 英文键（报告逻辑与语言解耦）
    KEYMAP = {"长度": "Length", "SS 类别": "SS class", "SS 置信": "SS conf",
              "TM 段": "TM", "拓扑": "Topology", "无序": "Disorder",
              "长程接触": "LR contact", "接触序": "Contact order",
              "ESM α": "ESM a", "ESM β": "ESM b", "基序": "Motifs",
              "同源等级": "Homology", "同源来源": "Source",
              "SP 命中": "SP hit", "SP 一致性": "SP ident",
              "SP 覆盖": "SP cov", "SP 蛋白名": "SP protein",
              "SP 基因": "SP gene", "SP 物种": "SP organism",
              "SP 证据": "SP PE", "事件种类": "Event type",
              "指标": "Metric"}
    for r in an:
        for z, e in KEYMAP.items():
            if z in r and e not in r:
                r[e] = r[z]
    for r in ng:
        for z, e in KEYMAP.items():
            if z in r and e not in r:
                r[e] = r[z]
    for r in scr:
        for z, e in KEYMAP.items():
            if z in r and e not in r:
                r[e] = r[z]

    doc = Document()
    st = doc.styles["Normal"]
    st.font.name = "Times New Roman"
    st.font.size = Pt(12)
    st.element.rPr.rFonts.set(qn("w:eastAsia"), "SimSun")

    def para(text, bold=False, size=12, align=None, after=6):
        p = doc.add_paragraph()
        r = p.add_run(text)
        r.bold = bold
        r.font.size = Pt(size)
        r.font.name = "Times New Roman"
        if align:
            p.alignment = align
        p.paragraph_format.space_after = Pt(after)
        return p

    def table(rows, size=7.5, max_rows=None):
        if max_rows:
            rows = rows[:max_rows]
        t = doc.add_table(rows=len(rows), cols=len(rows[0]))
        t.style = "Table Grid"
        for i, row in enumerate(rows):
            for j, cell in enumerate(row):
                c = t.cell(i, j)
                c.text = str(cell)
                for p in c.paragraphs:
                    for r in p.runs:
                        r.font.size = Pt(size)
                        r.font.name = "Times New Roman"
                        r.bold = (i == 0)
        doc.add_paragraph()

    def figure(name, width=6.2):
        f = os.path.join(out_dir, name)
        if os.path.exists(f):
            doc.add_picture(f, width=Inches(width))
            doc.paragraphs[-1].alignment = WD_ALIGN_PARAGRAPH.CENTER

    def pick(keys_en, keys_zh):
        """按语言取列名。"""
        for k in (keys_zh if zh else keys_en):
            if k in an[0] if an else False:
                return k
        return (keys_zh if zh else keys_en)[0]

    # ---------------- 统计量
    main_rows = [r for r in ng if str(r.get("Event type", "")).strip()
                 and not str(r.get("Event type", "")).startswith("  ")]
    col_total_pep = next((c for c in (ng[0].keys() if ng else [])
                          if c.endswith("Total pep") or "合计" in c and "肽" in c),
                         None)
    def total_pep(row):
        if col_total_pep:
            return _num(row.get(col_total_pep))
        vals = [_num(v) for k, v in row.items() if "肽" in k or "pep" in k.lower()]
        return max(vals) if vals else 0.0
    n_pep = sum(total_pep(r) for r in main_rows)
    inter = next((r for r in main_rows
                  if str(r.get("Event type", "")).startswith(("基因间区", "Intergenic"))), {})
    corr = next((r for r in main_rows
                 if str(r.get("Event type", "")).startswith(("已注释", "Annotated"))), {})
    def total_prot(row):
        keys = [k for k in row if "蛋白" in k or "prot" in k.lower()]
        vals = [_num(row.get(k)) for k in keys]
        return max(vals) if vals else 0.0
    # 同源分级列：英/中表头两种写法，缺一个就会全 0（曾因此把 6 个命中写成 0）
    def _hom(r):
        return (r.get("Homology class") or r.get("Homology")
                or r.get("同源分级") or "")
    hom = Counter(_hom(r) for r in an)
    sp_hits = [r for r in an if r.get("SP hit")]
    ref_strong = sum(1 for r in an if _hom(r) in ("strong", "identical"))

    # ---------------- 物种与数据规模（不写死本项目）
    _SPNAME = {"TP": "Treponema pallidum", "Mab": "M. abscessus",
               "Mch": "M. chelonae"}
    _codes = sorted({str(r.get("Species", "")).strip() for r in an
                     if str(r.get("Species", "")).strip()})
    _labels = [_SPNAME.get(c, c) for c in _codes] or ["the studied organism(s)"]
    if len(_labels) == 1:
        sp_title = _labels[0]
    elif len(_labels) == 2:
        sp_title = "%s and %s" % (_labels[0], _labels[1])
    else:
        sp_title = ", ".join(_labels[:-1]) + " and " + _labels[-1]
    n_raw = 0
    _ri = os.path.join(root, "00_input", "raw_info.tsv")
    if os.path.exists(_ri):
        with open(_ri, newline="", encoding="utf-8-sig") as fh:
            n_raw = max(sum(1 for _ in fh) - 1, 0)
    n_groups = len([c for c in (scr[0].keys() if scr else []) if c != "Metric"])

    # ---------------- 标题
    para((f"新基因事件鉴定报告：{sp_title} 镜像酶切蛋白质基因组学"
          if zh else
          f"Novel Gene Events in {sp_title}: Mirror-Digestion Proteogenomics"),
         bold=True, size=14, align=WD_ALIGN_PARAGRAPH.CENTER, after=4)
    para(f"ORFscape 自动生成｜{datetime.date.today():%Y-%m-%d}", size=10,
         align=WD_ALIGN_PARAGRAPH.CENTER, after=14)

    # ---------------- 摘要
    para("摘要" if zh else "Abstract", bold=True, size=12)
    if zh:
        para(f"本研究以 {n_raw} 个 raw（{n_groups} 个样品组）为输入，"
             f"经多层数据库建库、两段漏斗搜库、isoPSV 竞争仲裁与基因组坐标判定，"
             f"获得可信新肽段 {int(n_pep):,} 条；据此鉴定新编码基因 "
             f"{int(total_prot(inter)):,} 个（肽段 {int(total_pep(inter)):,} 条）"
             f"与已注释基因校正 {int(total_prot(corr)):,} 个（肽段 "
             f"{int(total_pep(corr)):,} 条）。对全部新编码蛋白完成：①多参数结构特征"
             f"（二级结构、跨膜拓扑、卷曲螺旋、低复杂度、无序度、ESM-2 注意力接触"
             f"统计）；②Swiss-Prot 与 RefSeq 同源比对及功能转移；③基序扫描。"
             f"Swiss-Prot 层面强同源 {hom.get('strong', 0) + hom.get('identical', 0)} 个、"
             f"中等 {hom.get('moderate', 0)} 个、弱 {hom.get('weak', 0)} 个、"
             f"无命中 {hom.get('none', 0)} 个，后者构成潜在全新编码序列候选。")
    else:
        para(f"From {n_raw} raw files ({n_groups} sample group(s)), a multi-layer "
             f"database, two-stage funnel searches, "
             f"isoPSV competitive arbitration and genome-coordinate assignment yielded "
             f"{int(n_pep):,} credible novel peptides, supporting {int(total_prot(inter)):,} "
             f"novel coding genes ({int(total_pep(inter)):,} peptides) and "
             f"{int(total_prot(corr)):,} annotated-gene corrections "
             f"({int(total_pep(corr)):,} peptides). All novel proteins were profiled for "
             f"structural features (secondary structure, TM topology, coiled-coil, "
             f"low-complexity, disorder, ESM-2 attention-contact statistics), annotated by "
             f"Swiss-Prot and RefSeq homology transfer, and screened for sequence motifs. "
             f"Against Swiss-Prot: {hom.get('strong', 0) + hom.get('identical', 0)} strong, "
             f"{hom.get('moderate', 0)} moderate, {hom.get('weak', 0)} weak and "
             f"{hom.get('none', 0)} without any hit, the latter being candidates for truly "
             f"novel coding sequences.")

    # ---------------- 方法
    para("1  数据与方法" if zh else "1  Data and methods", bold=True, size=12)
    if zh:
        para("1.1  事件分类：新肽段按六框 ORF 的坐标与链映射为基因组闭区间，仅与"
             "同链注释 CDS 比较；肽段与注释 CDS 同读框（相位校正后偏移可被 3 整除）"
             "判为 N 端延伸（已注释基因校正），其余（不同读框、反义链、基因间区）"
             "一律判为新编码基因；统计仅用通过过滤的可信新肽。", size=12)
        para("1.2  结构特征：理化性质由 Biopython ProtParam 计算；二级结构用简化"
             "Chou-Fasman（3 残基平滑）并给出置信度；跨膜段用 Kyte-Doolittle 19-aa "
             "滑窗（阈值 1.6），并按正电荷内侧规则给出 in/out 拓扑串；卷曲螺旋用七肽"
             "重复（a/d 位疏水）滑窗；低复杂度取 20-aa 窗口最小 Shannon 熵；无序度按"
             "TOP-IDP 倾向残基占比。", size=12)
        para("1.3  接触图：ESM-2（35M，本地权重）取末层注意力并按头平均，构建 "
             "L×L 接触倾向矩阵，统计长程接触密度、接触序以及接触基 α（i,i+3/4）与 "
             "β（i,i+2）分数。说明：本地权重不含监督的 contact-regression 头，"
             "故该矩阵为注意力倾向而非校准接触概率，仅作粗结构信号。", size=12)
        para("1.4  同源与功能：blastp（BLAST+ 2.15.0，E≤1e-3，最多 5 靶标）分别"
             "比对 Swiss-Prot（reviewed，本地快照）与本/交叉物种 RefSeq 蛋白组；"
             "按一致性与覆盖度分级（一致 ≥99%/≥99%、强 ≥70%/≥70%、中等 ≥30%/≥50%、"
             "弱、无），功能名由 Swiss-Prot 策展条目转移。另做 PROSITE 风格基序扫描"
             "（糖基化、P-loop、锌指、脂修饰、NLS 等，仅作提示）。", size=12)
    else:
        para("1.1  Event classification: each novel peptide was mapped through its "
             "six-frame ORF coordinates and strand to a genomic interval and compared "
             "only with same-strand annotated CDSs. Peptides sharing the annotated "
             "reading frame (phase-corrected offset divisible by 3) were called "
             "N-terminal extensions (annotation corrections); all others (different "
             "frame, antisense, intergenic) were called new coding genes. Only credible "
             "peptides passing filters were counted.", size=12)
        para("1.2  Structural features: physicochemical parameters from Biopython "
             "ProtParam; secondary structure by a simplified Chou-Fasman with 3-residue "
             "smoothing and a confidence score; TM segments by a Kyte-Doolittle 19-aa "
             "window (threshold 1.6) with an in/out topology string from the "
             "positive-inside rule; coiled-coil by heptad (a/d hydrophobicity) scanning; "
             "low complexity as the minimum 20-aa Shannon entropy; disorder as the "
             "TOP-IDP propensity fraction.", size=12)
        para("1.3  Contact maps: ESM-2 (35M, local weights) last-layer attention was "
             "averaged over heads to build an L×L contact-propensity matrix; long-range "
             "contact density, contact order and contact-based alpha (i,i+3/4) and beta "
             "(i,i+2) fractions were derived. Note: the local checkpoint lacks the "
             "supervised contact-regression head, so this matrix is an attention "
             "propensity rather than a calibrated contact probability and is used only "
             "as a coarse structural signal.", size=12)
        para("1.4  Homology and function: blastp (BLAST+ 2.15.0, E<=1e-3, max 5 targets) "
             "against Swiss-Prot (reviewed, local snapshot) and against the own and "
             "cross-species RefSeq proteomes; hits graded by identity/coverage "
             "(identical >=99/99, strong >=70/70, moderate >=30/50, weak, none). "
             "Function names were transferred from curated Swiss-Prot entries. PROSITE-"
             "style motif scanning (glycosylation, P-loop, zinc finger, lipid "
             "modification, NLS) is reported as hints only.", size=12)

    # ---------------- 结果
    para("2  新肽段筛选与新基因事件" if zh else
         "2  Peptide screening and novel gene events", bold=True, size=12)
    if scr:
        para("表 1  筛选漏斗（分样本，样本归属非互斥）" if zh else
             "Table 1  Screening funnel (per sample; non-exclusive attribution)",
             size=10)
        table([[r.get("Metric", "")] + [v for k, v in r.items() if k != "Metric"]
               for r in scr], size=7)
    if ng:
        para("表 2  事件分类（仅可信新肽）" if zh else
             "Table 2  Event classes (credible peptides only)", size=10)
        table([list(ng[0].values())] + [list(r.values()) for r in ng], size=6.8)
    figure("NG1_events_bar.png")
    figure("NG2_orf_distributions.png")

    # ---------------- 结构
    para("3  新编码蛋白的结构特征" if zh else
         "3  Structural features of novel proteins", bold=True, size=12)
    lens = [_num(r.get("Length")) for r in an]
    mws = [_num(r.get("MW(kDa)")) for r in an]
    pis = [_num(r.get("pI")) for r in an]
    grav = [_num(r.get("GRAVY")) for r in an]
    lrs = [_num(r.get("LR contact")) for r in an if r.get("LR contact")]
    tm_pos = sum(1 for r in an if _num(r.get("TM")) >= 1)
    ss_conf = [_num(r.get("SS conf")) for r in an if r.get("SS conf")]
    if zh:
        para(f"{len(an)} 个新编码蛋白：长度中位 {_median(lens):.0f} aa"
             f"（{min(lens):.0f}–{max(lens):.0f}），分子量中位 {_median(mws):.1f} kDa，"
             f"pI 中位 {_median(pis):.2f}，GRAVY 中位 {_median(grav):+.2f}；"
             f"{tm_pos} 个含至少 1 段预测跨膜螺旋；二级结构置信度中位 "
             f"{_median(ss_conf):.2f}。", size=12)
        if lrs:
            para(f"ESM-2 注意力接触统计：长程接触密度中位 {_median(lrs):.3f}"
                 f"（范围 {min(lrs):.3f}–{max(lrs):.3f}），接触序中位 "
                 f"{_median([_num(r.get('Contact order')) for r in an if r.get('Contact order')]):.2f}；"
                 f"接触基 α 分数中位 "
                 f"{_median([_num(r.get('ESM a')) for r in an if r.get('ESM a')]):.2f}、"
                 f"β 分数中位 "
                 f"{_median([_num(r.get('ESM b')) for r in an if r.get('ESM b')]):.2f}。",
                 size=12)
        ss = Counter(r.get("SS class", "") for r in an)
        para("二级结构类别：" + "；".join(f"{k} {v} 个" for k, v in ss.most_common())
             + "。", size=12)
    else:
        para(f"{len(an)} novel proteins: median length {_median(lens):.0f} aa "
             f"({min(lens):.0f}-{max(lens):.0f}), median MW {_median(mws):.1f} kDa, "
             f"median pI {_median(pis):.2f}, median GRAVY {_median(grav):+.2f}; "
             f"{tm_pos} carry at least one predicted TM helix; median SS confidence "
             f"{_median(ss_conf):.2f}.", size=12)
        if lrs:
            para(f"ESM-2 attention statistics: median long-range contact density "
                 f"{_median(lrs):.3f} ({min(lrs):.3f}-{max(lrs):.3f}); median contact "
                 f"order {_median([_num(r.get('Contact order')) for r in an if r.get('Contact order')]):.2f}; "
                 f"contact-based alpha median "
                 f"{_median([_num(r.get('ESM a')) for r in an if r.get('ESM a')]):.2f} and "
                 f"beta median "
                 f"{_median([_num(r.get('ESM b')) for r in an if r.get('ESM b')]):.2f}.",
                 size=12)
        ss = Counter(r.get("SS class", "") for r in an)
        para("Secondary-structure classes: "
             + "; ".join(f"{k} {v}" for k, v in ss.most_common()) + ".", size=12)
    figure("AN1_properties.png")
    figure("AN3_features.png")
    figure("AN4_contact_maps.png", width=6.0)

    # ---------------- 同源
    para("4  同源性比对与功能注释" if zh else
         "4  Homology and functional annotation", bold=True, size=12)
    if zh:
        para(f"Swiss-Prot（reviewed）层面：一致 {hom.get('identical', 0)}、"
             f"强 {hom.get('strong', 0)}、中等 {hom.get('moderate', 0)}、"
             f"弱 {hom.get('weak', 0)}、无命中 {hom.get('none', 0)}；"
             f"共 {len(sp_hits)} 个蛋白获得 Swiss-Prot 命中并有策展功能名。"
             f"RefSeq 本体层面另有 {ref_strong} 个强同源，多为旁系同源拷贝或注释"
             f"边界差异，其新基因身份需谨慎表述。", size=12)
    else:
        para(f"Against Swiss-Prot (reviewed): identical {hom.get('identical', 0)}, "
             f"strong {hom.get('strong', 0)}, moderate {hom.get('moderate', 0)}, "
             f"weak {hom.get('weak', 0)}, no hit {hom.get('none', 0)}; "
             f"{len(sp_hits)} proteins received a curated Swiss-Prot function. "
             f"At the RefSeq-self level a further {ref_strong} strong hits appear, mostly "
             f"paralogous copies or annotation-boundary variants whose novelty should be "
             f"stated cautiously.", size=12)
    if sp_hits:
        para("表 3  Swiss-Prot 命中与功能转移" if zh else
             "Table 3  Swiss-Prot hits and transferred function", size=10)
        heads = (["ORF", "物种", "长度", "一致性%", "覆盖%", "Swiss-Prot 蛋白名",
                  "基因", "物种来源", "PE"] if zh else
                 ["ORF", "Species", "Length", "Identity%", "Coverage%",
                  "Swiss-Prot protein", "Gene", "Organism", "PE"])
        rows = [heads] + [[r.get("ORF", ""), r.get("Species", ""),
                           r.get("Length", ""), r.get("SP ident", ""),
                           r.get("SP cov", ""), r.get("SP protein", "")[:40],
                           r.get("SP gene", ""), r.get("SP organism", "")[:28],
                           r.get("SP PE", "")] for r in sp_hits]
        table(rows, size=6.5, max_rows=40)
    figure("AN2_homology.png")
    # 基序
    motifs = Counter(m for r in an for m in str(r.get("Motifs", "")).split("; ")
                     if m)
    if motifs:
        para("表 4  基序扫描（提示性）" if zh else "Table 4  Motif hits (hints)",
             size=10)
        table([[("基序" if zh else "Motif"), ("蛋白数" if zh else "Proteins")]]
              + [[k, str(v)] for k, v in motifs.most_common(12)], size=7)

    # ---------------- 代表蛋白卡片
    para("5  代表性新编码蛋白" if zh else "5  Representative novel proteins",
         bold=True, size=12)
    cards = sorted([r for r in an if r.get("Homology class") == "none"],
                   key=lambda r: -_num(r.get("Length")))[:5]
    for r in cards:
        if zh:
            txt = (f"{r.get('ORF','')}（{r.get('Species','')}，{r.get('Length','')} aa）："
                   f"pI {r.get('pI','')}、GRAVY {r.get('GRAVY','')}；二级结构 "
                   f"{r.get('SS class','')}（α {r.get('α%','')}% / β {r.get('β%','')}%）；"
                   f"跨膜段 {r.get('TM','')}（拓扑 {r.get('Topology','') or '—'}）；"
                   f"无序 {r.get('Disorder','')}；长程接触 {r.get('LR contact','')}；"
                   f"无 Swiss-Prot/RefSeq 同源命中")
        else:
            txt = (f"{r.get('ORF','')} ({r.get('Species','')}, {r.get('Length','')} aa): "
                   f"pI {r.get('pI','')}, GRAVY {r.get('GRAVY','')}; SS class "
                   f"{r.get('SS class','')} (alpha {r.get('α%','')}% / beta "
                   f"{r.get('β%','')}%); TM {r.get('TM','')} (topology "
                   f"{r.get('Topology','') or '-'}); disorder {r.get('Disorder','')}; "
                   f"LR contact {r.get('LR contact','')}; no Swiss-Prot/RefSeq homolog")
        para(txt, size=11)

    # ---------------- 局限
    para("6  局限与建议" if zh else "6  Limitations and recommendations",
         bold=True, size=12)
    if zh:
        para("（1）N 端边界：新编码 ORF 的起始密码子来自六框预测，仅部分 ORF 有 "
             "N 端肽段或 N 端修饰证据；确证 5′ 边界需 N 端富集蛋白质组学"
             "（COFRADIC/TAILS）或 Ribo-seq 起始位点定位。", size=12)
        para("（2）结构预测：本报告未做三维结构建模；ESM-2 注意力接触为倾向性"
             "信号（无监督接触头），如需结构级结论应在具备 AlphaFold/ESMFold "
             "环境时补算并核对 pLDDT。", size=12)
        para("（3）功能注释：全部来自同源转移与基序提示，非实验验证；Swiss-Prot "
             "命中为本地快照，随时间推移可能新增条目。", size=12)
        para("（4）同源等级中的“强/一致”提示部分新 ORF 可能是旁系同源或注释边界"
             "差异，建议逐条核查其基因座上下文与转录证据后再作新基因主张。", size=12)
    else:
        para("(1) N-terminal boundaries: start codons come from six-frame prediction; "
             "only part of the ORFs have N-terminal peptides or N-terminal "
             "modifications. Confirming 5' boundaries needs N-terminomics "
             "(COFRADIC/TAILS) or Ribo-seq start-site mapping.", size=12)
        para("(2) Structure: no 3D modelling was performed; ESM-2 attention contacts are "
             "a propensity signal (no supervised contact head). Structure-level claims "
             "require AlphaFold/ESMFold runs with pLDDT inspection.", size=12)
        para("(3) Function: all annotations are homology/motif transfers, not "
             "experimental; the Swiss-Prot snapshot is local and grows over time.", size=12)
        para("(4) Strong/identical homology classes indicate some ORFs may be paralogs "
             "or annotation-boundary variants; verify locus context and transcript "
             "evidence before claiming a new gene.", size=12)

    para("7  产物文件" if zh else "7  Output files", bold=True, size=12)
    para("NG_table.tsv/xlsx、AN_table.tsv/xlsx（annotation + homologs）、"
         "NOV_screening_table.tsv、NG1/NG2 图、AN1–AN4 图、各系列结论句 md，"
         "以及 blast/ 下 Swiss-Prot 与 RefSeq 命中原始表。",
         size=11)

    name = "NovelGene_Report_zh.docx" if zh else "NovelGene_Report.docx"
    out = os.path.join(out_dir, name)
    doc.save(out)
    log(f"  报告已生成：{name}")
    return out
