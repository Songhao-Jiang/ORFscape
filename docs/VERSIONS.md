# Software versions

**English summary.** ORFscape itself is pure Python; everything that decides a
number (search engine, arbitration engine, predictors, gene prediction, splice
alignment, BLAST) is external and is *recorded, not assumed*. `orfscape doctor`
writes the exact versions it finds. Reports must carry the table below; a result
computed under a different major version of the arbitration engine or the search
engine is a different result.

## 1. ORFscape

| 项 | 版本 |
|---|---|
| ORFscape | **0.1.0**（论文发表前冻结在这条版本线；`python -m orfscape --version`） |
| Python | 3.9+，**实测 3.10** |
| 依赖（纯 Python 包） | `numpy`、`matplotlib`、`python-docx`、`Pillow`；GUI 仅 `tkinter`（随 Python 自带） |
| 平台 | Windows 10/11（完整功能，含 GUI 与 `rawdump`）；Python 包本身可移植 |

## 2. 决定数字的外部组件

| 组件 | 版本 | 用在哪 | 记录点 |
|---|---|---|---|
| **pFind** | **3.1**（含 pParse） | E1 / E2 / P6b 搜库；PSM FDR 1% | 项目 `validate_summary.md` 表头 + `env_report.tsv` |
| **isoPSV_multiengine** | 原核参考线 **2.0.177**；真核参考线 **2.0.184** | P7 竞争判读 | 由安装路径自动取版本写入汇总表头 |
| 谱图预测器 | AlphaPeptDeep（主）、MS2PIP、pDeep2、pDeep3 | isoPSV 四指标 | isoPSV 汇总表 |
| 格式转换 | msconvert（ProteoWizard）或 ThermoRawFileParser | raw → MGF | `env_report.tsv` |
| R | isoPSV 自带的 `ReadMGF_script.R` 所需（Rscript 在 PATH） | P7 抽谱 | `env_report.tsv` |
| **Prodigal** | 2.6.x（GPL-3，不随包分发） | 原核 L1 ab initio | `dbbuild` 日志；缺失时显式降级 |
| GeneMarkS-2 | 外部安装 | 原核 L1（可选第二源） | 同上 |
| **AUGUSTUS** | **3.5.0** + 物种模型（如 `cryptococcus_neoformans_neoformans_B`） | 真核第三期转录层 | `euk_stats.tsv` / 实验记录 |
| BRAKER3 | 仅生成脚本 | 真核注释（可选） | 不参与本版本数字 |
| **HISAT2 / StringTie** | 2.2.x（Linux，VM 内） | 真核 RNA-seq 组装 | `02_tx/<sp>/vm/` 日志 |
| orfipy | 外部安装 | L6 转录区 ORF | `L6_coords.tsv` 生成日志 |
| **BLAST+** | 实测 **2.11.0**（2.15 亦可） | L6 去冗余、新 ORF 同源 | 同源 README |
| Swiss-Prot | 2026-09 版 | `novannot` 同源分级 | `08_annotation/AN_table.tsv` 说明 |

## 3. 版本对数字的影响（必须交代的两处）

1. **isoPSV**：2.0.179–2.0.184 声明"零输出变化"，但 2.0.184 定版**必须**带
   `pfind_orphan_scope=global`，否则并行 pFind 探活按 result 目录匹配、cfg 却在
   param → 0/4 判死 EXIT1。2.0.177 及以前没有这个开关。
2. **pFind 3.1 + 半特异模板**：参数是口径的一部分（见 [CONSTRAINTS.md](CONSTRAINTS.md) §3.1），
   换搜索引擎（Comet/MSFragger/Sage）等于换结果，不可与既有 KEEP 数并列。

## 4. 怎么核对

```bash
python -m orfscape doctor --workdir orfscape_projects/<proj>
# → 04a_qc?/env_report.tsv：每个外部工具的名字、版本、状态（OK / missing）
```

外部工具的定位顺序：环境变量（`ORFSCAPE_PFIND` / `ORFSCAPE_ISOPSV` / `ORFSCAPE_BLAST` …）
→ `orfscape/paths.local.json`（git-ignored）→ 常见安装路径 → `PATH`。
仓库里没有任何绝对路径。

## 5. 写进报告的最小版本块

```
ORFscape 0.1.0 · pFind 3.1 · isoPSV_multiengine v2.0.184 · AUGUSTUS 3.5.0 ·
HISAT2/StringTie 2.2.x · BLAST+ 2.11.0 · Swiss-Prot 2026-09
```

（原核项目把 AUGUSTUS/HISAT2 换成 Prodigal / GeneMarkS-2。）
