# ORFscape pipeline — overview

**English summary.** ORFscape turns a genome into a layered ORF database and uses
DDA mass spectra to decide which ORFs were translated. Two locked reference lines
exist: **prokaryote P0→P9** (six layers, validated on *M. abscessus* / *M. chelonae*
and *T. pallidum*) and **eukaryote S1→S11** (spliced layers, validated on
*C. neoformans* H99). Automation is `orfscape run`: one child process per stage, one
automatic gate after each stage, a structured event stream the GUI can attach to,
and resume from any stage. Details: [PROKARYOTE.md](PROKARYOTE.md),
[EUKARYOTE.md](EUKARYOTE.md); non-negotiable rules: [CONSTRAINTS.md](CONSTRAINTS.md);
versions: [VERSIONS.md](VERSIONS.md).

本文是总入口。**两条线共用一个尾巴**（P6 四滤 → P6b 回搜 → P7 判读 → P8 事件/注释 →
P9 图/报告），差别只在"库怎么建、事件怎么归"。

---

## 1. 两条线一张图

```
            原核 P0–P9                                  真核 S1–S11
  raw ─ scan ─ genome ─ split ─┐          raw ─ scan ─ genome ─ split ─┐
                               │                                     │
        L0 注释 / L1 ab initio │              L0 剪接蛋白 / LJR 连接肽 │
        L2 六框 / L3 N 端      ├─ dbbuild     L2 转录本ORF(三期)       ├─ eukbuild
        L4/L5 / L6 转录证据    │              L1_SF 基因间六框(可选)   │  + dbbuild
                               │                                     │
                    E1 已注释库搜（40 任务）        E1 已注释库搜（40 任务）
                    QC / predbench                QC
                    E2 发现库搜（多库并行）         E2 发现库搜（6f/6g/8a/9i）
                               │                                     │
        ┌──────────────────────┴─────────────────────────────────────┐
        │  P6 四滤 → P6b 候选小库回搜 → P7 isoPSV 判读 → P8 事件/注释   │
        │  → P9 图/报告 → pack 结果包                                  │
        └────────────────────────────────────────────────────────────┘
```

| | 原核 | 真核 |
|---|---|---|
| 基因结构 | 连续 CDS | 断裂基因，必须按转录本拼接 |
| 层来源 | 基因组六框 / 上游起始 / 转录区 | 转录本三框 / 连接肽 / 从头预测 |
| 事件分类 | 四类互斥（N 端 > 不同框 > 反义 > 基因间） | 八类（另含内含子、C 端、同框内部、ncRNA） |
| 参考实测 | NTM 双物种：KEEP 2,421 行 / 1,471 肽；TP：KEEP 258 肽 / 402 行 | H99：KEEP 1,629 行 / 1,006 肽，事件 999 肽 / 537 ORF |
| 文档 | [PROKARYOTE.md](PROKARYOTE.md) | [EUKARYOTE.md](EUKARYOTE.md) |

## 2. 自动化契约（`orfscape run`）

```bash
python -m orfscape run --workdir orfscape_projects/<proj> --raw <raw 根> \
    --species-name "Treponema pallidum" [--from P6 --to P8c] [--skip P6b]
```

| 机制 | 行为 |
|---|---|
| 阶段表 | `orfscape/orchestrator.py` 的 `STEPS`：P0 → P1 → D0 → C → E1 → QC → E2 → P6 → P6b → P7 → P8a → P8b → P8c |
| 一阶段一子进程 | 子进程自己写事件流，父进程逐行转发日志 |
| **每阶段一道门** | 调用 `orfscape.checks`（`run_check`）；**FAIL 立即停**，不把空结果往下传 |
| 结构化事件 | `00_logs/progress.jsonl`：`stage.start / stage.end / stage.check`，GUI 只读这个文件 |
| 断点续跑 | `--from` / `--to` / `--skip`；已有产物按完成标记（glob）跳过 |
| 中断 | GUI 或人写 `00_logs/stop.flag`，编排器在行边界收尾 |
| 预估 | 上一轮实测阶段耗时写入 `00_logs/calibration.json`，下一轮 ETA 用真值 |
| 路径安全 | 项目根 `realpath` + 单段文件名白名单 + `dirname` 必须等于基目录 |

**已经自动的**：原核 13 个阶段全部；P6→P9 尾部两条线共用。
**还没进 `run` 的**（真核前缀段）：`eukbuild` / `euktx` / `sframe --build` 需要项目级参数
（基因组、GFF、转录本、preset），目前由脚本或人按 EUKARYOTE.md 逐条执行，再交给 `run --from E1`。
这是已知缺口，不假装已自动化。

### 每阶段的自检入口

```bash
python -m orfscape.checks --workdir <proj>              # 全阶段
python -m orfscape.checks --workdir <proj> --stage C    # 单阶段
```

`integrations/scripts/stage_check.py` 与 `preflight.py` 是给 agent 的薄壳，判据相同。

## 3. 阶段—产物—门的对照（原核）

| Stage | 命令 | 门（自动） |
|---|---|---|
| P0 | `scan` | raw > 0，scans > 0 |
| P1 | `genome` | fna/gff/faa 存在且非空 |
| D0 | `split` | MGF 数 == raw 数 |
| C | `dbbuild` | L0/L2 非空；有 L1 时必须有 `L1_coords.tsv` |
| E1 | `e1` | 任务数 == MGF 数；PSM > 0 |
| QC | `qc` / `qcfigs` / `e1figs` | 图表齐全 |
| PB | `predbench` | `04a_qc/predictor_summary.tsv` 有 `chosen` |
| E2 | `e2` | 任务数 == MGF × 库数；PSM > 0 |
| P6 | `novel` | 行数 > 0；层命中与各层 faa 非空一致 |
| P6b | `research` | 合并 ≥ P6；候选库 ≪ 全库 |
| P7 | `validate` | KEEP > 0；每条 MISSING 有原因 |
| P8 | `novgene` / `novannot` / `novreport` | 事件 > 0；无"未定位" |
| P9 | `novelfigs` / `pack` | 漏斗与参数图存在；`result_<YYYYMMDD>/` 结构完整 |

真核 S1–S11 的同表见 [EUKARYOTE.md](EUKARYOTE.md) §2。

## 4. 交付契约

每次任务开始日命名：`result_YYYYMMDD/`（英文）+ `result_YYYYMMDD_chinese/`（中文）。
包内 `report/ tables/ figures/ pdf/ README.md`；README 必须写清 KEEP 行/唯一肽、事件数、
注释数与口径。**PPT 不进结果包**，单独出中英两版。结果包用 `scripts/cneo_pack.py`
（先清空再拷）或项目自带的 packer；通用 `orfscape pack` **不覆盖已存在文件**，
复用旧目录会把上一轮的图与报告留在包里。

**版本号冻结**：论文正式发表前一律 **v0.1.0**——功能、修复、文档都进这一条线，不升号；
发布包名随之固定为 `ORFscape_v<版本>.zip`，GitHub 上只保留一个 **v0.1.0** Release/tag。
`build_release.py` 只打包 git 跟踪的文件，并拒绝个人路径、凭据与人名缩写；
`docs/VERSIONS.md` 记录与该版本一起冻结的外部组件版本。

## 5. 出问题时先看哪里

1. `00_logs/progress.jsonl` 最后一个 `stage.check` → 哪道门 FAIL；
2. 该阶段子进程日志尾部（编排器会打印）；
3. `python -m orfscape doctor --workdir <proj>` → 外部工具是否被找到、版本对不对；
4. `docs/CONSTRAINTS.md` 的"反模式"表 → 是否踩了已知静默失败。
