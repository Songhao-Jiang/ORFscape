# ORFscape v0.1 Release Notes

**v0.1.0 · 2026-11-13 · first public release**

ORFscape turns a genome into a layered ORF database and uses DDA mass spectra to
decide which ORFs were actually translated. This first release covers the
prokaryotic pipeline end to end, the eukaryotic layer builder, a desktop GUI, and
agent integration.

## Highlights

- **Six-layer database** — L0 annotated, L1 ab initio (Prodigal / GeneMarkS-2 /
  Augustus), L2 six-frame, L3 N-terminal extension, L4/L5 autolysis and
  contaminants, L6 transcript evidence, plus splicing junction peptides for
  eukaryotes. Every layer keeps its own header prefix so a novel peptide can be
  traced back to the layer that produced it.
- **P0→P9 orchestration with a gate after every stage** — `orfscape run` spawns one
  child process per stage, streams a structured event stream to
  `00_logs/progress.jsonl`, runs the stage checker, and stops on the first failed
  gate instead of carrying an empty result downstream. Resumable from any stage
  (`--from`/`--to`/`--skip`), interruptible through `00_logs/stop.flag`, and the
  next run's ETA is calibrated from the previous run's measured stage times.
- **Discovery search and arbitration** — E1 against the annotated proteome, E2
  against the discovery layers in parallel task trees, P6 novel-peptide filtering
  with an explicit layer-hit check, P6b candidate-ORF re-search against a small
  database, P7 per-spectrum arbitration with Tier1/Tier2/Tier3 and a card threshold
  recomputed from Tier1 each round.
- **Events and report** — mutually exclusive classification (N-terminal correction >
  different frame > antisense > intergenic), structure/function/homology annotation,
  Word report, funnel and parameter-comparison figures.
- **Desktop GUI** (tkinter, no third-party dependency) — English by default with a
  Chinese switch, light/dark themes, per-stage board with progress and ETA, run plan
  preview before starting, resumable failed stages, resource strip (CPU/memory/disk
  and process RSS), database panel with remote sizes, log search and right-click
  menu, export summary, taskbar progress (degrades to the title bar when the shell
  does not expose `ITaskbarList3`). The GUI only reads the event stream, so it can
  attach to a run started by the CLI or an agent, and closing it never stops the run.
- **Eukaryotic layer builder** — `orfscape eukbuild` splices CDS per transcript,
  reverse-complements minus-strand genes, reads each transcript's `transl_table`
  (1/2/3), emits spliced proteins, junction peptides (12 aa per side, `cut%3==0`,
  no stop) and an optional per-transcript three-frame ORF layer, then validates the
  result against the reference proteome.
- **Agent integration** — a shell-agnostic runbook plus adapters for ZCode, Claude
  Code and any generic chat model (`integrations/`).

## Validation in this release

| Organism | Proteins | Reference agreement | Junction peptides |
|---|---|---|---|
| *S. cerevisiae* R64 | 6,009 | 99.92% | 163 |
| *C. reinhardtii* (CC-503 v5.5) | 19,525 | 99.98% | 84,606 |
| *H. sapiens* (GRCh38.p14) | 145,213 | 98.23% exact / 99.22% with documented substitutions and non-ATG starts | 782,122 |

Human RefSeq `protein.faa` (136,807 records) and its I/L-deduplicated form
(90,373 unique) reproduce the published counts exactly; the differences from 100%
are N-terminal start codons (RefSeq forces M), genome-vs-RefSeq substitutions
recorded in the GFF, internal stops, and models without a RefSeq counterpart.
The test suite (7 entries, including a pixel-level GUI regression) passes on
Windows 10/11 with Python 3.10.

## Requirements

Python 3.9+ with `numpy`, `matplotlib`, `python-docx`; pFind 3.1 + pParse; msconvert
or ThermoRawFileParser; BLAST+; isoPSV; Prodigal (not bundled, GPL-3). Windows is
required for the bundled `rawdump`/`rawmethod` and the GUI's Win32 integrations.
Run `python -m orfscape doctor` to see what is present.

## Known limitations

- Eukaryotic **search** is not wired yet: the layers and the EUK gate exist, but
  E1/E2 assume the prokaryotic database layout.
- `LJR|` / `L2|` are not yet in the novel-peptide whitelist; novel peptides spanning
  an intron are not located by transcript coordinates.
- BRAKER3 support is script generation only.
- Card thresholds are recomputed per run; KEEP counts are not comparable across runs.

## License

MIT. Third-party components and their terms are listed in `NOTICE`.

---

# ORFscape v0.1 发布说明

**v0.1.0 · 2026-11-13 · 首次公开发布**

ORFscape 把基因组铺成多层 ORF 库，再用 DDA 质谱判断哪些 ORF 真的被翻译。本版包含
原核全流程、真核建库、桌面图形界面与 agent 接入。

## 主要内容

- **六层库**：L0 注释、L1 ab initio（Prodigal/GeneMarkS-2/Augustus）、L2 六框、
  L3 N 端延伸、L4/L5 自切与污染、L6 转录本证据；真核另加剪接连接肽层。每层保留
  自己的表头前缀，新肽可回溯到产生它的层。
- **P0→P9 编排，每阶段一道门**：`orfscape run` 逐阶段起子进程，事件流写
  `00_logs/progress.jsonl`，跑完立即自检；门不过即停，不把空结果带到下游。支持
  任意阶段续跑（`--from/--to/--skip`）、`stop.flag` 中断、按上一轮实测耗时标定 ETA。
- **发现搜库与仲裁**：E1 注释库、E2 发现库（并行任务树）、P6 新肽筛选（带层命中
  校验）、P6b 候选 ORF 小库回搜、P7 逐谱仲裁（Tier1/2/3，卡值每轮按 Tier1 重算）。
- **事件与报告**：互斥事件归类（N 端校正 > 不同框 > 反义 > 基因间区）、结构/功能/
  同源注释、Word 报告、漏斗与参数比较图。
- **桌面界面**（tkinter，无第三方依赖）：英文默认可切中文、明暗主题、阶段看板与
  ETA、开跑前的计划预览、失败阶段续跑、资源条（CPU/内存/磁盘与进程 RSS）、数据库
  面板与远端体积、日志搜索与右键菜单、导出汇总、任务栏进度（不可用时降级到标题栏）。
  界面只读事件流，因此能接管 CLI 或 agent 起的批次，关闭界面不影响跑批。
- **真核建库**：`orfscape eukbuild` 按转录本拼接 CDS、负链反向互补、逐转录本读
  `transl_table`（1/2/3），产出剪接蛋白层、连接肽层（两侧各 12 aa）、可选的转录本
  三框层，并与参考蛋白集做一致性校验。
- **Agent 接入**：与壳无关的 runbook，加 ZCode / Claude Code / 通用模型三套适配。

## 本版验证

| 物种 | 蛋白数 | 与参考一致 | 连接肽 |
|---|---|---|---|
| 酿酒酵母 R64 | 6,009 | 99.92% | 163 |
| 莱茵衣藻 CC-503 v5.5 | 19,525 | 99.98% | 84,606 |
| 人 GRCh38.p14 | 145,213 | 完全一致 98.23% / 计入参考替换与非 ATG 起始 99.22% | 782,122 |

人 RefSeq `protein.faa`（136,807 条）与 I/L 去重后的 90,373 条均与公开口径一致；
与 100% 的差异来自 N 端起始密码子（RefSeq 强制写 M）、GFF 记录的基因组与参考差异、
内部终止、以及无 RefSeq 对应的模型。7 个测试入口（含像素级界面回归）在 Windows 10/11
+ Python 3.10 全绿。

## 已知边界

- 真核**搜库**尚未接线：层与 EUK 门已就位，E1/E2 仍按原核库布局。
- `LJR|` / `L2|` 未进新肽层白名单；跨内含子新肽尚未按转录本坐标定位。
- BRAKER3 仅生成脚本。
- 卡值每轮重算，KEEP 数跨轮不可比。

## 许可

MIT。第三方组件与许可条款见 `NOTICE`。
