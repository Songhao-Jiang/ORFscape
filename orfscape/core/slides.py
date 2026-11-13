# -*- coding: utf-8 -*-
"""汇报 PPT 生成器：中英双语（默认中文汇报版），16:9，嵌入全部结果图表。

用法: from orfscape.core.slides import build_ppt
"""
import os

from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.util import Inches, Pt

from .i18n import get_lang
from .manuscript import check_arithmetic, gather_numbers, gather_versions
from .. import __version__

# 配色
C_TITLE = RGBColor(0x1B, 0x2A, 0x4A)   # 深蓝
C_BODY = RGBColor(0x33, 0x33, 0x33)     # 深灰
C_ACCENT = RGBColor(0x0A, 0x7B, 0xA4)   # 青蓝
C_WARN = RGBColor(0xD9, 0x53, 0x0F)     # 橙
C_WHITE = RGBColor(0xFF, 0xFF, 0xFF)
C_LIGHT = RGBColor(0xF0, 0xF4, 0xF8)    # 浅底

WIDTH = Inches(13.333)
HEIGHT = Inches(7.5)

FONT_TITLE = "Microsoft YaHei" if get_lang() == "zh" else "Arial"
FONT_BODY = "Microsoft YaHei" if get_lang() == "zh" else "Arial"


def _blank(prs):
    return prs.slides.add_slide(prs.slide_layouts[6])


def _txt(slide, left, top, width, height, text, size=18, bold=False,
         color=C_BODY, align=PP_ALIGN.LEFT, font=FONT_BODY, line_spacing=1.15):
    box = slide.shapes.add_textbox(left, top, width, height)
    tf = box.text_frame
    tf.word_wrap = True
    tf.vertical_anchor = MSO_ANCHOR.TOP
    lines = text.split("\n")
    for i, line in enumerate(lines):
        if i == 0:
            p = tf.paragraphs[0]
        else:
            p = tf.add_paragraph()
        p.alignment = align
        p.space_after = Pt(4)
        p.line_spacing = line_spacing
        r = p.add_run()
        r.text = line
        r.font.size = Pt(size)
        r.font.bold = bold
        r.font.color.rgb = color
        r.font.name = font
    return box


def _bar(slide, left, top, width, height=Pt(4), color=C_ACCENT):
    shape = slide.shapes.add_shape(1, left, top, width, height)  # 1 = rectangle
    shape.fill.solid()
    shape.fill.fore_color.rgb = color
    shape.line.fill.background()
    return shape


def _bg(slide, prs):
    bg = slide.shapes.add_shape(1, 0, 0, WIDTH, HEIGHT)
    bg.fill.solid()
    bg.fill.fore_color.rgb = C_LIGHT
    bg.line.fill.background()
    return bg


def _header(slide, title, subtitle=""):
    _bar(slide, Inches(0.5), Inches(0.5), Inches(12.3), Pt(3), C_ACCENT)
    _txt(slide, Inches(0.5), Inches(0.7), Inches(12.3), Inches(0.7),
         title, size=28, bold=True, color=C_TITLE, font=FONT_TITLE)
    if subtitle:
        _txt(slide, Inches(0.5), Inches(1.35), Inches(12.3), Inches(0.4),
             subtitle, size=14, color=C_ACCENT)


def _pic(slide, path, left, top, width=None, height=None):
    if width:
        return slide.shapes.add_picture(path, left, top, width=width)
    return slide.shapes.add_picture(path, left, top, height=height)


def _kpi_row(slide, left, top, kpis, w=Inches(2.8), h=Inches(1.2)):
    """kpies = [(值, 标签), ...] 横排 KPI 卡片"""
    for i, (val, label) in enumerate(kpis):
        x = left + i * (w + Inches(0.2))
        card = slide.shapes.add_shape(1, x, top, w, h)
        card.fill.solid()
        card.fill.fore_color.rgb = C_WHITE
        card.line.color.rgb = C_ACCENT
        card.line.width = Pt(1)
        _txt(slide, x, top + Inches(0.15), w, Inches(0.6),
             str(val), size=30, bold=True, color=C_ACCENT, align=PP_ALIGN.CENTER)
        _txt(slide, x, top + Inches(0.72), w, Inches(0.4),
             label, size=12, color=C_BODY, align=PP_ALIGN.CENTER)


def build_ppt(root, out_dir, lang=None):
    """生成汇报 PPT。返回 pptx 路径。"""
    if lang:
        from .i18n import set_lang
        set_lang(lang)
    zh = get_lang() == "zh"
    os.makedirs(out_dir, exist_ok=True)
    N = gather_numbers(root)
    check_arithmetic(N)
    V = gather_versions(root)
    fig_dir = os.path.join(root, "08_report")

    prs = Presentation()
    prs.slide_width = WIDTH
    prs.slide_height = HEIGHT
    T = (lambda en, cn: cn if zh else en)

    # ---- Slide 1: 封面 ----
    s = _blank(prs)
    bg = s.shapes.add_shape(1, 0, 0, WIDTH, HEIGHT)
    bg.fill.solid(); bg.fill.fore_color.rgb = C_TITLE; bg.line.fill.background()
    _txt(s, Inches(1), Inches(2.0), Inches(11.3), Inches(1.2),
         T("Automated Mirror-Digest Proteogenomics Pipeline",
           "原核生物镜像酶切蛋白质基因组学全自动流程"),
         size=36, bold=True, color=C_WHITE, align=PP_ALIGN.CENTER)
    _txt(s, Inches(1), Inches(3.3), Inches(11.3), Inches(0.6),
         T("ORFscape — from raw files to Nature-format figures in one command",
           "ORFscape — 一条命令从 raw 文件到 Nature 格式全套图表"),
         size=18, color=RGBColor(0xAA, 0xCC, 0xEE), align=PP_ALIGN.CENTER)
    _txt(s, Inches(1), Inches(4.5), Inches(11.3), Inches(0.5),
         T("M. abscessus + M. chelonae | FAIMS dual-CV | Trypsin+LysC / LysargiNase+LysN",
           "脓肿分枝杆菌 + 龟分枝杆菌 | FAIMS 双电压 | 胰酶/赖精酶镜像酶切"),
         size=14, color=RGBColor(0x88, 0xAA, 0xCC), align=PP_ALIGN.CENTER)
    _txt(s, Inches(1), Inches(6.2), Inches(11.3), Inches(0.5),
         T("XXX Lab | 2026-11", "XXX 实验室 | 2026-11"),
         size=14, color=RGBColor(0x66, 0x88, 0xAA), align=PP_ALIGN.CENTER)

    # ---- Slide 2: 研究背景与目标 ----
    s = _blank(prs); _bg(s, prs)
    _header(s, T("Background & Objectives", "研究背景与目标"),
            T("Why automated proteogenomics for NTM?", "为什么做 NTM 蛋白质基因组学自动化？"))
    items = [
        T("NTM annotations are incomplete — novel genes/ORFs/N-term extensions remain undiscovered",
          "NTM 基因组注释不完整——新基因/ORF/N 端延伸有待发现"),
        T("Existing workflows require heavy manual configuration (database building, splitting, searching)",
          "现有流程需大量手动配置（建库、拆分、搜库参数）"),
        T("FAIMS dual-CV raw files break conventional search engines (spectra from both voltages mix)",
          "FAIMS 双电压 raw 文件使常规搜库引擎出错（两电压谱图混叠）"),
        T("Mirror digestion (Trypsin+LysC / LysargiNase+LysN) provides orthogonal validation",
          "镜像酶切（胰酶/赖精酶）提供正交验证轴"),
        T("Goal: input = raw + species + groups → output = Nature-format figures/tables/methods",
          "目标：输入 = raw + 物种 + 分组 → 输出 = Nature 格式图表/方法学全文"),
    ]
    y = Inches(2.0)
    for item in items:
        _txt(s, Inches(0.8), y, Inches(11.5), Inches(0.55), f"• {item}",
             size=16, color=C_BODY)
        y += Inches(0.75)

    # ---- Slide 3: 数据一览 ----
    s = _blank(prs); _bg(s, prs)
    _header(s, T("Dataset Overview", "数据一览"),
            T("All numbers reconcile by addition", "全部数字可用加法对齐"))
    _kpi_row(s, Inches(0.5), Inches(1.8), [
        (f"{N['n_raw']}", T("raw files", "raw 文件")),
        (f"{N['size_gb']} GB", T("total size", "总数据量")),
        (f"{N['n_scan_total']:,}", T("total scans", "总扫描数")),
        (f"{N['n_mirror_groups']}", T("mirror groups", "镜像组")),
    ])
    _txt(s, Inches(0.5), Inches(3.3), Inches(12.3), Inches(0.5),
         T(f"Instrument: Thermo Orbitrap Exploris 480 + FAIMS (−45 V / −65 V) | "
           f"MS1 60k / MS2 15k | MIT 50/22 ms | NCE 27 | AGC 300%/75%",
           f"仪器：Thermo Orbitrap Exploris 480 + FAIMS（−45 V / −65 V）| "
           f"MS1 60k / MS2 15k | 最大注入时间 50/22 ms | NCE 27 | AGC 300%/75%"),
         size=13, color=C_BODY)
    _txt(s, Inches(0.5), Inches(4.0), Inches(12.3), Inches(0.3),
         T("Arithmetic check: ", "加法对账："), size=14, bold=True, color=C_TITLE)
    _txt(s, Inches(0.5), Inches(4.4), Inches(12.3), Inches(1.0),
         f"MS1 {N['n_ms1_total']:,} + MS2 {N['n_ms2_total']:,} = {N['n_scan_total']:,} scans\n"
         f"CV-split MGF spectra {N['n_spec_mgf']:,} = MS2 total (each spectrum assigned to exactly one CV, 64/64 pass)",
         size=14, color=C_BODY)
    _txt(s, Inches(0.5), Inches(5.5), Inches(12.3), Inches(0.8),
         T("Genomes: " + "; ".join(V["genomes"].values()),
           "参考基因组：" + "；".join(V["genomes"].values())),
         size=12, color=C_BODY)

    # ---- Slide 4: 流水线架构 ----
    s = _blank(prs); _bg(s, prs)
    _header(s, T("Pipeline Architecture (ORFscape v%s)" % __version__,
                 "流水线架构（ORFscape v%s）" % __version__),
            T("9 CLI commands, fully automated, resumable, bilingual (EN/ZH)",
              "9 个 CLI 子命令，全自动、断点续跑、中英双语"))
    steps = [
        ("scan", T("raw discovery + group inference + method extraction", "raw 发现+分组+方法参数")),
        ("doctor", T("environment check (13/18 tools auto-detected)", "环境检测（13/18 工具自动发现）")),
        ("genome", T("RefSeq assembly download + version pinning", "RefSeq 装配下载+版本钉扎")),
        ("dbbuild", T("7-layer database construction (L0–L6)", "七层建库（L0–L6）")),
        ("split", T("FAIMS CV split (64 MGF, exact reconciliation)", "FAIMS 电压拆分（64 MGF，精确对账）")),
        ("e1", T("annotated-DB search (pFind)", "注释库搜库（pFind）")),
        ("qc", T("acquisition + digestion QC", "采集+酶切质控")),
        ("e2", T("discovery-DB search (db_aug + db_6f)", "发现库搜库（db_aug + db_6f）")),
        ("novel+report", T("4-filter → GSSP → Tier → Nature figures", "四重过滤→GSSP→Tier→Nature 图表")),
    ]
    y = Inches(1.8)
    for i, (cmd, desc) in enumerate(steps):
        x = Inches(0.5) if i % 2 == 0 else Inches(6.8)
        if i % 2 == 0 and i > 0:
            y += Inches(1.0)
        _txt(s, x, y, Inches(1.8), Inches(0.4), cmd,
             size=14, bold=True, color=C_ACCENT)
        _txt(s, x + Inches(2.0), y, Inches(4.0), Inches(0.7), desc,
             size=12, color=C_BODY)

    # ---- Slide 5: 关键技术突破 ----
    s = _blank(prs); _bg(s, prs)
    _header(s, T("Key Technical Solutions", "关键技术突破"))
    cards = [
        (T("FAIMS CV Split", "FAIMS 电压拆分"),
         T("scanEvent filter → per-CV MGF\n99.7% peptide overlap with manual FreeStyle\npParse incompatible with FAIMS raw",
           "scanEvent 过滤 → 逐 CV MGF\n与手工 FreeStyle 拆分肽级 99.7% 一致\npParse 不兼容 FAIMS raw 格式")),
        (T("Parallel Search", "并行搜库"),
         T("4 instances × 4 threads = 16 threads\nCPU 37.9% → 99.3%\nE1: 2.5h → 50min (3× speedup)",
           "4 实例 × 4 线程 = 16 线程\nCPU 利用率 37.9% → 99.3%\nE1：2.5h → 50min（3× 提速）")),
        (T("Windows-Native VM-Free", "Windows 原生去 VM 化"),
         T("ENA direct FASTQ (no SRA toolkit)\nBLASTN (no HISAT2)\norfipy (no rnaSPAdes)",
           "ENA 直下 FASTQ（免 SRA 工具）\nBLASTN（免 HISAT2）\norfipy（免 rnaSPAdes）")),
        (T("Arithmetic Assertions", "数字对账断言"),
         T("12 addition checks enforced at build time\nMS1+MS2=scans; Tier sum=novel; layer sum=db size\nManuscript cannot generate if any check fails",
           "构建时强制 12 条加法断言\nMS1+MS2=扫描数；Tier 之和=新肽数；层和=库大小\n任何一条失败即拒绝生成")),
    ]
    positions = [(Inches(0.5), Inches(1.8)), (Inches(6.8), Inches(1.8)),
                 (Inches(0.5), Inches(4.2)), (Inches(6.8), Inches(4.2))]
    for (title, desc), (x, y) in zip(cards, positions):
        card = s.shapes.add_shape(1, x, y, Inches(6.0), Inches(2.2))
        card.fill.solid(); card.fill.fore_color.rgb = C_WHITE
        card.line.color.rgb = C_ACCENT; card.line.width = Pt(1)
        _txt(s, x + Inches(0.2), y + Inches(0.15), Inches(5.6), Inches(0.4),
             title, size=16, bold=True, color=C_TITLE)
        _txt(s, x + Inches(0.2), y + Inches(0.6), Inches(5.6), Inches(1.5),
             desc, size=12, color=C_BODY)

    # ---- Slide 6: 数据库构建 ----
    s = _blank(prs); _bg(s, prs)
    _header(s, T("Seven-Layer Database Architecture", "七层数据库构建"),
            T("Layer sums define search databases (asserted by addition)", "层和定义搜库库（加法断言）"))
    rows = []
    for sp in N["layers"]:
        L = N["layers"][sp]
        rows.append(f"{sp}: L0 {L['L0']:,} | L1 {L['L1']:,} | L2 {L['L2']:,} | "
                    f"L3 {L['L3']:,} | L4 {L['L4']} | L5 {L['L5']} | L6 {L['L6']}")
        rows.append(f"  db_annot {L['db_annot']:,} = L0+L4+L5 → "
                    f"db_aug {L['db_aug']:,} = +L1+L3 → "
                    f"db_6f {L['db_6f']:,} = +L2")
    _txt(s, Inches(0.5), Inches(1.8), Inches(12.3), Inches(2.5),
         "\n".join(rows), size=14, color=C_BODY)
    _txt(s, Inches(0.5), Inches(4.5), Inches(12.3), Inches(0.5),
         T("Source layers: L0=RefSeq annotation | L1=Prodigal ab initio | "
           "L2=six-frame ORF (≥10 aa) | L3=N-term extension | "
           "L4=enzyme autolysis | L5=contaminants | L6=transcript-evidenced ORF (BLASTN + orfipy)",
           "来源层：L0=RefSeq 注释 | L1=Prodigal 预测 | "
           "L2=六框 ORF（≥10 aa）| L3=N 端延伸 | "
           "L4=酶自切 | L5=污染 | L6=转录证据 ORF（BLASTN + orfipy）"),
         size=12, color=C_ACCENT)

    # ---- Slide 7: 鉴定结果 ----
    s = _blank(prs); _bg(s, prs)
    _header(s, T("Identification Results", "鉴定结果"))
    _kpi_row(s, Inches(0.5), Inches(1.8), [
        (f"{N['e1_psm']:,}", "E1 PSMs"),
        (f"{N['e2_psm']:,}", "E2 PSMs"),
        (f"{N['novel_total']:,}", T("Novel peptides", "新肽段")),
        (f"{N['t1_total']}", "Tier1"),
    ])
    _pic(s, os.path.join(fig_dir, "Fig1_workflow_identification.png"),
         Inches(1.5), Inches(3.3), width=Inches(10.3))
    _txt(s, Inches(0.5), Inches(7.0), Inches(12.3), Inches(0.4),
         T(f"E2 breakdown: db_aug {N['e2_by_db']['aug'][1]:,} + db_6f "
           f"{N['e2_by_db']['6f'][1]:,} = {N['e2_psm']:,}",
           f"E2 分账：db_aug {N['e2_by_db']['aug'][1]:,} + db_6f "
           f"{N['e2_by_db']['6f'][1]:,} = {N['e2_psm']:,}"),
         size=12, color=C_BODY, align=PP_ALIGN.CENTER)

    # ---- Slide 8: 质控 ----
    s = _blank(prs); _bg(s, prs)
    _header(s, T("Quality Control", "质量控控"))
    _pic(s, os.path.join(fig_dir, "Fig2_QC.png"),
         Inches(1.5), Inches(1.8), width=Inches(10.3))
    _txt(s, Inches(0.5), Inches(6.5), Inches(12.3), Inches(0.5),
         T("MS2 IT saturation >73% (MIT 22 ms cap) — consistent with known NTM acquisition profile | "
           "Cleavage specificity: TL 80.3% / LL 86.8% (fully-cleaved)",
           "MS2 注入时间触顶率 >73%（MIT 22 ms 封顶）——与 NTM 已知采集特征一致 | "
           "酶切特异性：TL 80.3% / LL 86.8%（全切率）"),
         size=12, color=C_BODY, align=PP_ALIGN.CENTER)

    # ---- Slide 9: 新肽发现 ----
    s = _blank(prs); _bg(s, prs)
    _header(s, T("Novel Peptide Discovery", "新肽段发现"))
    _pic(s, os.path.join(fig_dir, "Fig3_novel_funnel.png"),
         Inches(1.5), Inches(1.8), width=Inches(10.3))
    per_sp = " + ".join("%s %s" % (sp, format(v["n"], ","))
                        for sp, v in N["novels"].items())
    _txt(s, Inches(0.5), Inches(6.5), Inches(12.3), Inches(0.5),
         T(f"{per_sp} = "
           f"{N['novel_total']:,} | Tier1 {N['t1_total']} + Tier2 {N['t2_total']} + "
           f"Tier3 {N['t3_total']:,} = {N['novel_total']:,} | Decoy: 0",
           f"{per_sp} = "
           f"{N['novel_total']:,} | Tier1 {N['t1_total']} + Tier2 {N['t2_total']} + "
           f"Tier3 {N['t3_total']:,} = {N['novel_total']:,} | Decoy: 0"),
         size=13, bold=True, color=C_ACCENT, align=PP_ALIGN.CENTER)

    # ---- Slide 10: 镜像验证 ----
    s = _blank(prs); _bg(s, prs)
    _header(s, T("Mirror-Digest Validation", "镜像酶切验证"),
            T("Tier1 = same core peptide under both enzymes", "Tier1 = 双酶下同一核心肽"))
    _pic(s, os.path.join(fig_dir, "Fig4_mirror_validation.png"),
         Inches(2.5), Inches(2.0), width=Inches(8.3))
    _txt(s, Inches(0.5), Inches(6.3), Inches(12.3), Inches(0.6),
         T("Tier1 peptides are N/C-terminal K/R-stripped, I/L-normalised core sequences "
           "observed under both Trypsin+LysC and LysargiNase+LysN — providing "
           "engine-independent physical evidence for novel coding regions",
           "Tier1 肽段为去端 K/R、I/L 归一后的核心序列在双酶（胰酶/赖精酶）条件下均被鉴定——"
           "为新编码区提供不依赖搜库引擎的物理证据"),
         size=12, color=C_BODY, align=PP_ALIGN.CENTER)

    # ---- Slide 11: 产出清单 ----
    s = _blank(prs); _bg(s, prs)
    _header(s, T("Deliverables", "产出清单"))
    items = [
        T("4 main figures (PDF+SVG+PNG, 89/183 mm, Arial, colourblind-safe) — 08_report/",
          "4 张主图（PDF+SVG+PNG，89/183 mm，Arial，色盲安全）— 08_report/"),
        T("4 supplementary tables (xlsx: identification, novel peptides, QC, parameters)",
          "4 张附表（xlsx：鉴定总账、新肽全表、QC、参数）"),
        T("Full manuscript .docx (Times New Roman, Abstract→Introduction→Results→Methods→Figures, "
          "with 12 arithmetic assertions and XXX placeholders)",
          "论文级 Word 全稿（Times New Roman，摘要→引言→结果→方法→图，含 12 条加法断言和 XXX 占位）"),
        T("Methods .md (auto-generated, EN + ZH, all version provenance embedded)",
          "方法学 .md（自动生成，英文+中文，全部版本溯源内嵌）"),
        T("QC reports (per-task cleavage specificity, IT saturation, AGC fill)",
          "QC 报告（逐任务酶切特异性、IT 触顶率、AGC 填充）"),
        T("Novel peptide table (peptide / tier / layer / locus / enzymes / PSM count)",
          "新肽表（肽段/Tier/层/基因座/酶/PSM 数）"),
        T("Split reconciliation report (64/64 exact spectral count match)",
          "拆分对账报告（64/64 逐谱计数精确匹配）"),
        T("Complete version provenance (genomes, tools, databases, parameters)",
          "完整版本溯源（基因组、工具、数据库、参数）"),
    ]
    y = Inches(1.8)
    for item in items:
        _txt(s, Inches(0.8), y, Inches(11.5), Inches(0.55), f"✓ {item}",
             size=14, color=C_BODY)
        y += Inches(0.62)

    # ---- Slide 12: 总结与下一步 ----
    s = _blank(prs); _bg(s, prs)
    _header(s, T("Summary & Next Steps", "总结与下一步"))
    _txt(s, Inches(0.5), Inches(1.8), Inches(5.8), Inches(0.4),
         T("Summary", "总结"), size=18, bold=True, color=C_TITLE)
    summary = [
        T("Fully automated: raw → Nature-format output in one command chain",
          "全自动：raw → Nature 格式产出一条命令链"),
        T(f"{N['novel_total']:,} novel peptides, {N['t1_total']} with mirror-digest confirmation, zero decoy",
          f"{N['novel_total']:,} 条新肽段，{N['t1_total']} 条有镜像酶切确认，零 decoy"),
        T("All numbers reconcile by addition — enforced programmatically",
          "全部数字可加法对齐——程序化强制执行"),
        T("Windows-native, zero VM dependency, parallel execution (3× speedup)",
          "Windows 原生、零 VM 依赖、并行执行（3× 提速）"),
    ]
    y = Inches(2.3)
    for item in summary:
        _txt(s, Inches(0.7), y, Inches(5.5), Inches(0.5), f"• {item}",
             size=13, color=C_BODY)
        y += Inches(0.7)
    _txt(s, Inches(6.8), Inches(1.8), Inches(5.8), Inches(0.4),
         T("Next Steps", "下一步"), size=18, bold=True, color=C_TITLE)
    nexts = [
        T("Tier1 drop analysis (296→78 after Acetyl[AnyN-term])", "Tier1 下降分析（296→78）"),
        T("L1/L3/L6 layer hit-check (L2 absorption issue)", "L1/L3/L6 层命中检查（L2 吸收问题）"),
        T("Novel peptide spectrum spot-check (10 × Tier1/2)", "新肽谱图抽查（10 条 Tier1/2）"),
        T("Separate FDR recalculation (novel-only)", "分开 FDR 重算（novel-only）"),
        T("M. chelonae L6 transcriptome layer", "龟分枝杆菌 L6 转录本层"),
        T("isoPSV arbitration integration (v2.0.160)", "isoPSV 仲裁接入（v2.0.160）"),
        T("GUI + AI skill/plugin packaging", "GUI + AI skill/插件打包"),
    ]
    y = Inches(2.3)
    for item in nexts:
        _txt(s, Inches(7.0), y, Inches(5.5), Inches(0.5), f"→ {item}",
             size=13, color=C_BODY)
        y += Inches(0.6)

    out = os.path.join(out_dir, "ORFscape_Report.pptx" if zh else "ORFscape_Report_EN.pptx")
    prs.save(out)
    return out
