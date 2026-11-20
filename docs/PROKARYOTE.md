# Prokaryote line — P0 → P9 (locked)

**English summary.** Six layers from one genome, two searches (annotated, discovery),
four filters, a candidate-only re-search, per-spectrum arbitration, four mutually
exclusive event classes, annotation, figures and a result pack. Every stage has a
gate; `orfscape run` chains them. Reference runs: NTM (*M. abscessus* + *M. chelonae*,
KEEP 2,421 rows / 1,471 peptides) and TP (*T. pallidum*, KEEP 258 peptides / 402 rows).

命令里的 `<P>` = `orfscape_projects/<proj>`。逐步命令与门见 §2，口径见
[CONSTRAINTS.md](CONSTRAINTS.md)，版本见 [VERSIONS.md](VERSIONS.md)。

## 1. 层与库

| 层 | 前缀 | 内容 | 来源 |
|---|---|---|---|
| L0 | `L0_` | 注释蛋白 | GFF/GenBank CDS（按 `transl_table` 翻译） |
| L1 | `L1_` | ab initio 预测 | Prodigal / GeneMarkS-2 的 `.faa` 放进 `02_db/<sp>/` |
| L2 | `L2_` | 六框 ORF（stop-to-stop） | 基因组 |
| L3 | `L3\|` | N 端上游同框起始 | 注释 + 基因组 |
| L4 / L5 | `L4_` / `L5_` | 自切肽 / 常见污染 | 内置 |
| L6 | `L6\|` | 转录区内的 ORF | RNA-seq（见 §3） |

库：`db_annot` = L0+L4+L5（E1 用）；`db_aug` = L0+L1+L2+L3+L6（E2 用）；`db_6f{,1..3}`
= 大库按长度三池（E2 并行）。每个库配等量 `REV_` 诱饵。

> **前缀即溯源**：新肽落在哪一层，由 FASTA 头前缀决定。`L1_` / `L2_` / `L3|` / `L6|`
> 四个前缀必须同时进白名单，漏一个就是整层静默丢失。

## 2. 阶段表

| Stage | 命令 | 产物 | 门 |
|---|---|---|---|
| **P0 数据确认** | `python -m orfscape scan --workdir <P> --raw <raw 根> --species-name "…"` | `00_input/*` | raw > 0、scans > 0 |
| **P1 基因组** | `python -m orfscape genome --workdir <P> [--accession GCF_…] [--group bacteria]` | `01_genome/<sp>/` | fna/gff/faa 非空 |
| **P1b 参考库**（可选） | `python -m orfscape uniprot --workdir <P> --sp <sp>` | `uniprot_ref.faa` | 蛋白数 > 0 |
| **D0 拆分** | `python -m orfscape split --workdir <P>` | `03_mgf/*.mgf` | MGF 数 == raw 数 |
| **C 建库** | `python -m orfscape dbbuild --workdir <P>` | `02_db/<sp>/L0..L6 + db_*` | L0/L2 非空；有 L1 必有 `L1_coords.tsv` |
| **E1 已注释搜库** | `python -m orfscape e1 --workdir <P> --thread <T>` | `04_search_annot/` | 任务数 == MGF 数、PSM > 0 |
| **QC** | `python -m orfscape qc --workdir <P>`（+ `qcfigs` / `e1figs`） | `04a_qc/` | 图与表齐全 |
| **PB 预测器基线**（建议） | `python -m orfscape predbench --workdir <P> --n 600` | `04a_qc/predictor_*` | `chosen` 已定 |
| **E2 发现搜库** | `python -m orfscape e2 --workdir <P> --dbs aug,6f --workers <W> --thread <T>` | `04_search_disc/` | 任务数 == MGF × 库数、PSM > 0 |
| **P6 四滤** | `python -m orfscape novel --workdir <P> [--species <sp>]` | `05_novel/<sp>/novel_peptides.tsv` | 行数 > 0；层命中与各层 faa 一致 |
| **P6b 回搜** | `python -m orfscape research --workdir <P> --thread <T> --workers <W> --db-mode novel` | `05b_research/` | 合并 ≥ P6；候选库 ≪ 全库 |
| **P7 判读** | `python -m orfscape validate --workdir <P> --from-research [--cards s,c] [--reuse <dir>]` | `06b_validate_research/novel_validated.tsv` | KEEP > 0；MISSING 有 `v_miss_reason` |
| **P8a 事件** | `python -m orfscape novgene --workdir <P> --from 06b_validate_research` | `07_novgene/NG_table.tsv` | 事件 > 0、无"未定位" |
| **P8b 注释** | `python -m orfscape novannot --workdir <P> --from 06b_validate_research [--no-esm]` | `08_annotation/` | 属性/同源/特征表齐全 |
| **P8c 报告** | `python -m orfscape novreport --workdir <P> --from 06b_validate_research` | `08_report/*.docx` | 报告含非零统计 |
| **P9 图** | `python -m orfscape novelfigs --workdir <P> --from 06b_validate_research` | `08_report/` | 漏斗 + 参数比较图 |
| **交付** | `python scripts/<proj>_p13_*.py` + `pack` | `P13_results/`、`result_YYYYMMDD/` | 见 [PIPELINE.md](PIPELINE.md) §4 |

一步到位（可断点续跑）：

```bash
python -m orfscape run --workdir <P> --raw <raw 根> --species-name "<sp>"
```

## 3. L6 转录证据（可选层，必须在 C 之前）

配方（ENA 公开数据，实测 5.5 h → 16 min）：

1. ENA 里选 run：同株、读长合适、rRNA 占比低；
2. **块内去重后 BLAST**（同 block 75% 重复，先合并再去冗余）；
3. GFF 先把 rRNA 区掩膜；
4. 取 `±depth ≥ 3` 的连续转录区；
5. `orfipy` 出 BED，按基因组 `transl_table`（细菌 11）翻译；
6. CDS 与注释重叠 > 80% 的丢弃（那是已注释基因）；
7. 落盘 `02_db/<sp>/L6_tx_orfs.faa` + `L6_coords.tsv`，`dbbuild` 才会并进 `db_aug`。

```bash
python -m orfscape l6map --workdir <P> --sp <sp> --bed <orfipy.bed> --table 11
```

## 4. 两条线的共用尾巴

P6 四滤（消化集 / E1 实测集 / ≥7 aa 连续子串 / 等重变体）→ P6b 只放候选 ORF 的小库回搜
→ P7 isoPSV 竞争判读（Tier1/Tier2/Tier3 + 卡值）→ P8 事件与注释 → P9 图与报告。
细节与判据见 [CONSTRAINTS.md](CONSTRAINTS.md) §2，真核差异见 [EUKARYOTE.md](EUKARYOTE.md)。

## 5. 参考实测（用于对表，不作为阈值来源）

### NTM 双物种（*M. abscessus* + *M. chelonae*，11-17 定稿）

| 项 | 值 |
|---|---|
| 层 | L0 + L1（Prodigal + GeneMarkS-2）+ L2 六框 + L3 + L6（Mab 25,273 / Mch 12,752 ORF） |
| E1 / E2 | 全库搜，PSM 见项目报告 |
| P6 → P6b | 四滤后进回搜，候选库缩小两个数量级 |
| **P7 KEEP** | **2,421 行 / 1,471 唯一肽**（Mab 1,182 / Mch 1,239 行） |
| Tier | 974 / 1,039 / 408 |
| MISSING | 0 |
| 事件 | 1,039 肽 / 484 ORF（N 端另有 1,382 肽 / 249 ORF） |
| 注释 | 480 个 ORF |

### TP（*T. pallidum*，v6 定稿）

| 项 | 值 |
|---|---|
| 层 | L1 = Prodigal 61 + GeneMarkS-2 39；L6 已入库（43 ORF / 110 未注释转录区） |
| P6 / P6b | 823 行 / 670 肽 → 850 行 |
| **P7 KEEP** | **258 肽 / 402 行**；MISSING 4；unmapped 0 |
| 事件 | 111 ORF：N 端 25 / 基因间 22 / 不同框 23 / 反义 41 |
| 卡值 | cos 0.86（比 06b 更严仍净增 28 肽——L3 层以前没提交的真实回收） |

> **两轮卡值不可比**：NTM 与 TP 的 score/cos 各自按 Tier1 重新标定。报 KEEP 数时必须同时报卡值。

## 6. 已知坑（原核侧）

| 现象 | 真因 | 处理 |
|---|---|---|
| L6 整层 0 命中 | 白名单写成 `L6_`（实际是 `L6\|`） | 四个前缀一起检查 |
| 事件全归到反义 | bedtools 链向列取错 | 用 `-s` 显式指定，校验链向 |
| P6b 秒回"0 命中"、任务全 rc=1 | 注释蛋白 + 候选 ORF 同库 → pFind `Group.exe` 组合爆炸 | `--db-mode novel` |
| 换库后命中还是旧的 | 结果目录被复用 | `db.sha256` 指纹，不符挪 `.stale` 重搜 |
| N 端事件虚高 | 一条肽同读框多 CDS 重名被重复计数 | 按位点去重（同 contig/链/读框、重叠 > 80% 计一次） |
| `makeblastdb` 建库失败 | 目标盘不允许建目录 | 换本地盘（如 `D:/blastdb_*`） |
