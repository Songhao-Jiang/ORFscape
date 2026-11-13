"""sjh_ppt.py -- build the 0.SJH_results briefing deck (16:9, python-pptx).

Usage:
    python sjh_ppt.py --iso-img A.png --igv-img B.png --out C.pptx
"""
import argparse

from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.enum.shapes import MSO_SHAPE

NAVY = RGBColor(0x1F, 0x38, 0x64)
BLUE = RGBColor(0x2F, 0x6F, 0xBA)
GREY = RGBColor(0x59, 0x59, 0x59)
LIGHT = RGBColor(0xEA, 0xF0, 0xF8)
WHITE = RGBColor(0xFF, 0xFF, 0xFF)
RED = RGBColor(0xC0, 0x50, 0x4D)
ORANGE = RGBColor(0xE8, 0xA3, 0x3D)

W, H = 13.333, 7.5


def blank(prs):
    return prs.slides.add_slide(prs.slide_layouts[6])


def box(slide, x, y, w, h, fill=None, line=None, shape=MSO_SHAPE.RECTANGLE):
    sh = slide.shapes.add_shape(shape, Inches(x), Inches(y), Inches(w), Inches(h))
    if fill is None:
        sh.fill.background()
    else:
        sh.fill.solid()
        sh.fill.fore_color.rgb = fill
    if line is None:
        sh.line.fill.background()
    else:
        sh.line.color.rgb = line
        sh.line.width = Pt(1)
    sh.shadow.inherit = False
    return sh


def text(slide, x, y, w, h, s, size=18, color=GREY, bold=False, align=PP_ALIGN.LEFT,
         anchor=MSO_ANCHOR.TOP, spacing=1.15):
    tb = slide.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    tf = tb.text_frame
    tf.word_wrap = True
    tf.vertical_anchor = anchor
    tf.margin_left = tf.margin_right = 0
    tf.margin_top = tf.margin_bottom = 0
    lines = s.split("\n")
    for i, ln in enumerate(lines):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.alignment = align
        p.line_spacing = spacing
        r = p.add_run()
        r.text = ln
        r.font.size = Pt(size)
        r.font.bold = bold
        r.font.color.rgb = color
        r.font.name = "Microsoft YaHei"
    return tb


def header(slide, title, kicker=""):
    box(slide, 0, 0, W, 0.95, fill=NAVY)
    box(slide, 0, 0.95, W, 0.06, fill=ORANGE)
    text(slide, 0.55, 0.13, 11.5, 0.55, title, size=25, color=WHITE, bold=True)
    if kicker:
        text(slide, 0.55, 0.68, 11.5, 0.26, kicker, size=11.5,
             color=RGBColor(0xC9, 0xD6, 0xEA))


def footer(slide, n):
    text(slide, 0.55, H - 0.45, 6, 0.3,
         "ORFscape · 梅毒螺旋体未注释肽段 · 2026-11-09", size=9,
         color=RGBColor(0x9A, 0x9A, 0x9A))
    text(slide, W - 1.1, H - 0.45, 0.6, 0.3, str(n), size=10,
         color=RGBColor(0x9A, 0x9A, 0x9A), align=PP_ALIGN.RIGHT)


def kpi(slide, x, y, value, label, color=NAVY):
    box(slide, x, y, 2.55, 1.35, fill=LIGHT)
    text(slide, x, y + 0.18, 2.55, 0.6, value, size=30, color=color,
         bold=True, align=PP_ALIGN.CENTER)
    text(slide, x, y + 0.85, 2.55, 0.4, label, size=11, color=GREY,
         align=PP_ALIGN.CENTER)


def bullets(slide, x, y, w, items, size=15, gap=0.42):
    for i, (head, body) in enumerate(items):
        yy = y + i * gap * 1.55
        box(slide, x, yy + 0.09, 0.11, 0.11, fill=BLUE, shape=MSO_SHAPE.OVAL)
        text(slide, x + 0.28, yy, w - 0.28, 0.32, head, size=size, color=NAVY,
             bold=True)
        if body:
            text(slide, x + 0.28, yy + 0.36, w - 0.28, 0.5, body, size=11.5,
                 color=GREY)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--iso-img", required=True)
    ap.add_argument("--igv-img", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    prs = Presentation()
    prs.slide_width, prs.slide_height = Inches(W), Inches(H)
    n = 0

    # ---- 1 cover
    s = blank(prs)
    box(s, 0, 0, W, H, fill=NAVY)
    box(s, 0, 4.55, W, 0.07, fill=ORANGE)
    text(s, 1.0, 2.15, 11.3, 1.0, "梅毒螺旋体未注释肽段", size=40, color=WHITE,
         bold=True)
    text(s, 1.0, 3.15, 11.3, 0.8, "肽段质控 · 双库注释 · 新基因事件", size=24,
         color=RGBColor(0xC9, 0xD6, 0xEA))
    text(s, 1.0, 5.0, 11.3, 0.4, "Data: 498 novel peptides | ORFscape v0.1.0",
         size=14, color=RGBColor(0xC9, 0xD6, 0xEA))
    n += 1

    # ---- 2 inputs
    s = blank(prs)
    header(s, "输入：498 条未注释肽 + 两条 pFind 谱表", "单物种 · 镜像双酶 · 已有人工清单")
    kpi(s, 0.55, 1.35, "498", "新肽（T 291 / L 207）")
    kpi(s, 3.25, 1.35, "1,282", "匹配 PSM（T 777 / L 505）")
    kpi(s, 5.95, 1.35, "261", "涉及 ORF")
    kpi(s, 8.65, 1.35, "100", "带修饰肽")
    text(s, 0.55, 3.05, 12.2, 0.4, "参考与输入文件", size=15, color=NAVY, bold=True)
    rows = [
        "基因组：GCA_000410535.2（CP004010.2，1,139,633 bp，978 CDS）",
        "注释：genomic.gtf（GTF 2.2）→ 978 条 CDS，含 locus_tag / gene / product",
        "蛋白库：NCBI 978 条、UniProt UP000000811 1,027 条",
        "谱表：Proteome/Trypsin 57,627 PSM、Proteome/R29 28,393 PSM",
        "输入清单：梅毒重注释/Novel_pep.txt（Sequence / Modification / Proteins / Cleavage）",
    ]
    for i, r in enumerate(rows):
        text(s, 0.55, 3.5 + i * 0.52, 12.2, 0.4, "·  " + r, size=12.5, color=GREY)
    footer(s, n + 1)
    n += 1

    # ---- 3 pipeline
    s = blank(prs)
    header(s, "四步流程：先定位、再分级、后判读", "与 NTM 主线同一套新肽筛选口径")
    steps = [("① 基因组定位", "ORF 名 → 片段 → 负链反向互补\n→ 翻译 → I/L 归一串定位", BLUE),
             ("② Tier 分级", "Tier1 跨酶同核心\nTier2 同 ORF 双酶\nTier3 单酶", NAVY),
             ("③ isoPSV 收口", "全部肽过 isoPSV\n镜像证据可独立救回\nTier1 锚定卡值", ORANGE),
             ("④ 注释与出图", "双库四分类 → 事件归类\nGTF/BED + 两本 PDF", RED)]
    for i, (h1, h2, col) in enumerate(steps):
        x = 0.55 + i * 3.15
        box(s, x, 1.55, 2.85, 0.62, fill=col)
        text(s, x, 1.68, 2.85, 0.4, h1, size=15, color=WHITE, bold=True,
             align=PP_ALIGN.CENTER)
        box(s, x, 2.17, 2.85, 1.55, fill=LIGHT)
        text(s, x + 0.14, 2.33, 2.57, 1.3, h2, size=11.5, color=GREY)
        if i < 3:
            text(s, x + 2.85, 1.62, 0.3, 0.4, "→", size=20, color=ORANGE)
    text(s, 0.55, 4.05, 12.2, 0.4, "本批数据的两条特殊口径", size=15, color=NAVY,
         bold=True)
    bullets(s, 0.55, 4.55, 6.0, [
        ("没有一对肽在两酶中序列完全相同",
         "镜像酶切必然留不同末端 → Tier1 只能按核心区（去 K/R + I/L 归一）判定"),
    ])
    bullets(s, 6.9, 4.55, 5.9, [
        ("isoPSV 未过 ≠ 错误鉴定",
         "镜像证据是独立物理证据轴；Tier1/2 可由镜像单独保留（mirror_only）"),
    ])
    footer(s, n + 1)
    n += 1

    # ---- 4 mapping + tier
    s = blank(prs)
    header(s, "定位与分级：498/498 定位成功", "抽样 20 条用 ORF 全长翻译回验，20/20 一致")
    kpi(s, 0.55, 1.3, "498/498", "基因组定位成功")
    kpi(s, 3.2, 1.3, "149", "Tier1 跨酶同核心", RGBColor(0x2E, 0x75, 0xB6))
    kpi(s, 5.85, 1.3, "121", "Tier2 同 ORF 双酶", RGBColor(0x1F, 0x38, 0x64))
    kpi(s, 8.5, 1.3, "228", "Tier3 单酶", RGBColor(0x7F, 0x7F, 0x7F))
    text(s, 0.55, 3.05, 12.2, 0.4, "坐标换算", size=15, color=NAVY, bold=True)
    text(s, 0.55, 3.5, 6.0, 0.4, "正链：nt_start = orf_start + off × 3；nt_end = 起点 + 3L − 1",
         size=13, color=GREY)
    text(s, 0.55, 3.95, 6.0, 0.4, "负链：hi = orf_end − off × 3；nt_start = hi − 3L + 1",
         size=13, color=GREY)
    box(s, 6.95, 3.35, 5.8, 1.55, fill=LIGHT)
    text(s, 7.15, 3.5, 5.4, 1.3,
         "Tier1 用核心区的理由\n本批 498 条肽中，跨酶序列完全一致的对数为 0；\n"
         "去掉两端 K/R、I/L 归一后才有 149 对可比。",
         size=11.5, color=GREY)
    text(s, 0.55, 5.35, 12.2, 0.4, "ORF 命名：orf|frame|strand|start-end|，坐标 1-based 闭区间",
         size=12, color=RGBColor(0x8A, 0x8A, 0x8A))
    footer(s, n + 1)
    n += 1

    # ---- 5 db annotation
    s = blank(prs)
    header(s, "双库注释：UniProt 零增量，NCBI 是超集", "四分类口径（用户指定）")
    tab = [("db_class", "肽段", "ORF"),
           ("None（两库都没有）", "271", "213"),
           ("NCBI annotated", "222", "47"),
           ("UniProt annotated", "0", "0"),
           ("NCBI & UniProt annotated", "5", "5")]
    for i, (c1, c2, c3) in enumerate(tab):
        y = 1.45 + i * 0.62
        fill = NAVY if i == 0 else (WHITE if i % 2 else LIGHT)
        col = WHITE if i == 0 else GREY
        box(s, 0.55, y, 7.6, 0.6, fill=fill, line=RGBColor(0xD6, 0xDF, 0xEC))
        text(s, 0.75, y + 0.14, 4.6, 0.4, c1, size=13, color=col,
             bold=(i == 0))
        text(s, 5.4, y + 0.14, 1.2, 0.4, c2, size=13, color=col,
             bold=(i == 0), align=PP_ALIGN.RIGHT)
        text(s, 6.8, y + 0.14, 1.2, 0.4, c3, size=13, color=col,
             bold=(i == 0), align=PP_ALIGN.RIGHT)
    bullets(s, 8.55, 1.55, 4.3, [
        ("UniProt 零增量", "凡 UniProt 命中的，NCBI 都已收录"),
        ("5 条两库都有", "全部 7–8 aa 短肽，落在已知蛋白里的偶然短串"),
        ("建议", "按 db_class = 两库都有 且长度 ≤8 剔除 → 498 → 493"),
    ], size=14, gap=0.55)
    footer(s, n + 1)
    n += 1

    # ---- 6 isoPSV
    s = blank(prs)
    header(s, "isoPSV 收口：KEEP 281 / DROP 217", "判读 255/498 肽；卡值 score ≥0.32、cos ≥0.90")
    kpi(s, 0.55, 1.3, "255", "判读肽数（51%）")
    kpi(s, 3.2, 1.3, "281", "KEEP", RGBColor(0x2E, 0x75, 0xB6))
    kpi(s, 5.85, 1.3, "217", "DROP", RGBColor(0x9E, 0x4A, 0x45))
    tab2 = [("保留依据", "肽数"),
            ("isoPSV+mirror", "20"),
            ("mirror_only（isoPSV 未过，镜像保留）", "50"),
            ("isoPSV（Tier3 过三门）", "11"),
            ("no_isoPSV_mirror（无结论，Tier1/2 按镜像）", "200")]
    for i, (c1, c2) in enumerate(tab2):
        y = 3.05 + i * 0.52
        fill = NAVY if i == 0 else (WHITE if i % 2 else LIGHT)
        col = WHITE if i == 0 else GREY
        box(s, 0.55, y, 7.6, 0.5, fill=fill, line=RGBColor(0xD6, 0xDF, 0xEC))
        text(s, 0.75, y + 0.12, 6.0, 0.3, c1, size=12, color=col, bold=(i == 0))
        text(s, 6.9, y + 0.12, 1.1, 0.3, c2, size=12, color=col, bold=(i == 0),
             align=PP_ALIGN.RIGHT)
    bullets(s, 8.55, 3.15, 4.3, [
        ("mirror_only 不是缺陷", "48% 的保留来自镜像证据单独成立"),
        ("Tier3 三门通过率 8.3%", "与 NTM 侧的 8.4% 同量级"),
    ], size=13, gap=0.5)
    footer(s, n + 1)
    n += 1

    # ---- 7 events
    s = blank(prs)
    header(s, "新基因事件：268 条，落在 261 个 ORF", "只与同链 CDS 比较，相位 %3 判定是否同读框")
    kpi(s, 0.55, 1.3, "205", "基因间区新编码基因", RGBColor(0x2F, 0x6F, 0xBA))
    kpi(s, 3.2, 1.3, "63", "已注释基因 N 端延伸", RGBColor(0xC0, 0x50, 0x4D))
    kpi(s, 5.85, 1.3, "258", "间区事件涉及肽数")
    kpi(s, 8.5, 1.3, "240", "N 端延伸涉及肽数")
    text(s, 0.55, 3.0, 12.2, 0.4, "间区 205 条的细分", size=15, color=NAVY, bold=True)
    box(s, 0.55, 3.45, 5.9, 1.5, fill=LIGHT)
    text(s, 0.75, 3.6, 5.5, 1.3, "纯基因间区 125 条（164 肽）\n同链不同读框 80 条（94 肽）",
         size=13, color=GREY)
    box(s, 6.85, 3.45, 5.9, 1.5, fill=LIGHT)
    text(s, 7.05, 3.6, 5.5, 1.3,
         "63 条 N 端延伸对应 240 肽\n（多条肽落在同一基因的延伸区）",
         size=13, color=GREY)
    text(s, 0.55, 5.25, 12.2, 0.4,
         "事件级 268 条 vs 肽级 498 条：一个 ORF 的多条肽归并为一条事件。",
         size=12, color=RGBColor(0x8A, 0x8A, 0x8A))
    footer(s, n + 1)
    n += 1

    # ---- 8 mirror spectrum image
    s = blank(prs)
    header(s, "逐肽证据：镜像谱（实际 ↑ / 理论 ↓）", "isoPSV_spectra.pdf · 498 页 · 一条肽一页")
    s.shapes.add_picture(a.iso_img, Inches(0.5), Inches(1.25), width=Inches(7.6))
    bullets(s, 8.35, 1.45, 4.5, [
        ("四行表头", "谱名 / 肽段 / Judgment / 相似性"),
        ("离子梯", "b 下 y 上，切口带上下短横"),
        ("实测 vs 理论", "等高镜像，共用 0 基线"),
        ("ppm 图", "同横轴对齐，±20 ppm"),
    ], size=13, gap=0.62)
    footer(s, n + 1)
    n += 1

    # ---- 9 IGV image
    s = blank(prs)
    header(s, "位点示意：IGV 风格三轨图", "IGV_like_annotation.pdf · 268 页 · 一条事件一页")
    s.shapes.add_picture(a.igv_img, Inches(0.5), Inches(1.4), width=Inches(8.4))
    bullets(s, 9.2, 1.45, 3.7, [
        ("CDS 轨", "重叠基因自动分行\n标签两行：locus+gene / product"),
        ("新 ORF 轨", "红=N 端延伸、蓝=纯间区\n黄=同链异读框"),
        ("肽段轨", "矩形，宽度=基因组跨度"),
        ("窗口", "ORF ±25%（最少 ±1.5 kb）"),
    ], size=12.5, gap=0.85)
    footer(s, n + 1)
    n += 1

    # ---- 10 conclusions
    s = blank(prs)
    header(s, "结论", "")
    bullets(s, 0.75, 1.5, 11.8, [
        ("498 条肽全部定位成功，Tier1 149 / Tier2 121 / Tier3 228",
         "Tier1 只能用核心区判定——本批没有跨酶完全一致的肽对，这是镜像酶切的必然"),
        ("268 条新基因事件：间区 205 + N 端延伸 63，落在 261 个 ORF",
         "N 端延伸的 63 条事件覆盖 240 条肽，说明多个肽落在同一基因的延伸区"),
        ("UniProt 对本批零增量，NCBI 是超集",
         "双库注释的价值在本批体现为「NCBI 已注释 222 肽」，而非双库互补"),
        ("isoPSV 判读 51%，KEEP 281 · mirror_only 占 48%",
         "镜像证据是独立证据轴；把 mirror_only 当缺陷会误读这套设计"),
    ], size=15, gap=0.62)
    footer(s, n + 1)
    n += 1

    # ---- 11 limits
    s = blank(prs)
    header(s, "局限与下一步", "")
    bullets(s, 0.75, 1.45, 5.9, [
        ("isoPSV 覆盖面 51%", "450 条谱未产出结果；需查 1.unlabeled_result（短肽？NCE 参数？）"),
        ("5 条 7–8 aa 短肽疑似假新", "建议剔除后重算总数 498 → 493"),
    ], size=14, gap=0.62)
    bullets(s, 6.9, 1.45, 5.9, [
        ("与人工 xlsx 基线对照", "69 对镜像肽 ↔ Tier1 149；34 个新 ORF ↔ 268 事件"),
        ("29 条同源需 BLASTP", "本包未做"),
    ], size=14, gap=0.62)
    box(s, 0.55, 4.55, 12.2, 1.55, fill=LIGHT)
    text(s, 0.8, 4.72, 11.7, 1.3,
         "交付物：Table1/2/3 · GTF + BED · isoPSV_spectra.pdf（498 页） · "
         "IGV_like_annotation.pdf（268 页） · 本报告 md/docx · 实验记录 docx · 本 PPT",
         size=13, color=NAVY)
    footer(s, n + 1)

    prs.save(a.out)
    print("saved", a.out, "slides", len(prs.slides))


if __name__ == "__main__":
    main()
