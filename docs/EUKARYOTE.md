# Eukaryote line — S1 → S11 (locked)

**English summary.** Eukaryotic genes are spliced, so the layers change: spliced
protein (L0), junction peptides (LJR), transcript ORFs in three independent phases
(vendor assembly / VM re-assembly / AUGUSTUS ab initio), and an optional
intergenic six-frame layer gated by coding density (R = 6/d). The shared tail
(P6→P9) runs unchanged, but events are classified into eight classes at peptide
level with gene geometry taken from the **reference** annotation. Reference run:
*C. neoformans* H99 — E1 404,069 PSM, KEEP 1,629 rows / 1,006 peptides,
999 peptides / 537 ORFs, 22.3% species-specific at the amino-acid level.

真核线第一次接真实质谱是 H99（2026-11）。本文件是它的 as-built 记录，命令可直接搬到
另一个真菌/真核基因组，只换物种码与基因组三件套。

## 1. 与原核的四处硬差异

| 项 | 原核 | 真核 |
|---|---|---|
| 蛋白层 | 连续 CDS 翻译 | **按转录本拼接外显子**，密码子表读 GFF 的 `transl_table` |
| 连接肽 | 无 | **LJR**：每个剪接点两侧各 12 aa，`cut % 3 == 0` |
| 转录层 | 无（L6 是转录区 ORF） | **L2 三期互不混池**：厂家组装 / VM 组装 / AUGUSTUS 从头 |
| 事件 | 四类互斥 | **八类**，且基因几何**只能取参考注释** |

## 2. 阶段表（S1–S11）

`<P>` = `proj_<物种>` 项目根；所有 `--out` 必须是**绝对路径**（相对路径会落进 ORFscape 仓库）。

| Stage | 命令 | 产物 | 门 |
|---|---|---|---|
| **S1 数据入库** | `scan` → `genome --group fungi` → `split` | `00_input/`、`01_genome/`、`03_mgf/` | MGF 数 == raw 数；MS2 总数与 scan 表一致 |
| **S2 转录层（三期）** | ① `euktx --gtf vendor/merged.gtf` ② VM 里 `run_rna.sh`（HISAT2→StringTie→`--merge`）后 `euktx --gtf vm/merged.gtf` ③ `run_augustus.sh` + `cneo_l1_augustus.py` | `02_tx/<sp>/tx_model.gff`、`transcripts*.fa` | 三期条数与预期一致；AUGUSTUS 模型必须带 `CDS` 行与 `protein_coding` |
| **S3 建库** | 先过六框门 `sframe --sp <sp>`；需要时 `sframe --min-aa 60 --build`；然后 `eukbuild` → `dbbuild` | `db_annot`、`db_6f{,1..3}`、`db_6g{,1..3}`、`db_8a{,1..3}`、`db_9i{,1..3}` | L0 与参考蛋白逐条一致率 ≥ 99.5%；各池条数在预期 ±1% |
| **S4 E1** | `e1 --thread 16` → `qc` | `04_search_annot/` | 任务 40/40；PSM > 0；覆盖度写进 QC |
| **S5 E2** | `e2 --dbs 6f1..3` / `6g1..3` / `8a1..3` / `9i1..3` | `04_search_disc/` | 任务数 == MGF × 库数；PSM > 0 |
| **S6 P6 四滤** | `novel --workdir <P>` | `05_novel/<sp>/novel_peptides.tsv` | 行数 > 0，decoy 0；四滤的注释排除集**非空** |
| **S7 P6b** | `research --db-mode novel --workers 6 --thread 3` | `05b_research/` | 候选库缩小 ≥ 10×；合并 ≥ P6 |
| **S8 P7** | `validate --from-research --cards 0.31,0.9`（增量加 `--reuse`） | `06c_validate_full/` | KEEP > 0；MISSING 0 或有明确 `v_miss_reason`；**换肽集必删 `unlabeled/`** |
| **S9 P8** | `novgene --events euk`；`novannot --from …`；`cneo_report.py` | `07_novgene/euk_events.tsv`、`08_annotation/` | 八类事件齐；`intergenic` 不为 0；无未定位 |
| **S10 交付** | `cneo_cases.py` / `cneo_pdfs.py` / `cneo_report_ppt.py` / `cneo_pack.py` | `docs/pdf/`、`result_YYYYMMDD{,_chinese}/` | 页数 == 唯一肽/ORF 数；包内文件数固定；英文包 0 中文字符 |
| **S11 同源** | `cneo_homology.py {queries,blast,matrix,figs}` | `09_novel_homology/` | 特异性计数与交付表逐字节一致；`matrix` 可复算 |

## 3. 六框适用性门（建库前，先算再做）

真核基因组通常比原核稀疏，六框是否值得做只取决于一个数：

```
d = ΣCDS / 基因组                 # 编码密度
R = 六框残基 / 注释蛋白残基 ≈ 6 / d   # 库被放大多少倍
```

| 判定 | R | d | 做法 |
|---|---|---|---|
| compact | ≤ 15 | ≥ 40% | 可做（库与转录本层同量级） |
| moderate | 15–30 | 20–40% | 只在没有转录组时做 |
| sparse | 30–60 | 10–20% | 不建议 |
| extreme | > 60 | < 10% | 不做，改转录本 + 从头预测 |

另有算力闸：六框残基 > 150 M aa 或 ORF > 200 万条时不做。
判完再看**基因间占比**：整库 39 万里只有 8.9% 在注释基因外 → 只建基因间六框。

```bash
python -m orfscape sframe --workdir <P> --sp <sp>                    # 只判定，约 10 s
python -m orfscape sframe --workdir <P> --sp <sp> --min-aa 60 --build # 建基因间六框库
```

H99 实测：d = 65.0%、R = 9.2× → 可做；基因间 34,807/390,196（8.9%）→ 建 8,000 ORF。
E2 只出 316 PSM，P7 留 12 行且全部与组装转录本重叠 → **零新增基因**，这一层对 H99 是负结果。

## 4. 八类事件与基因几何

| 类别 | 判定 |
|---|---|
| intergenic | 肽不在任何**参考注释基因**跨度内 |
| intronic | 在基因内、CDS 外、不延伸至两端 |
| nterm_ext / cterm_ext | CDS 外且越过首/末端 |
| inframe_internal | CDS 内同框（ORF 有延伸、肽无） |
| alt_frame | CDS 范围内读框不同，或 ORF 延伸而肽无延伸 |
| antisense | 落在注释基因反义链 |
| ncRNA_translation | 落在注释 ncRNA 内 |

两条必须记住的规则：

1. **基因几何只取 `01_genome/<sp>/genomic.gff`**。用 `tx_model.gff` 会把 MSTRG/AUGUSTUS
   转录本按区间挂到注释基因上，把 span 撑进基因间区（`intergenic` 恒 0），无 CDS 的新模型
   还会被判成"非编码基因"（肽全记成 `ncRNA_translation`）。坐标映射仍用 `tx_model.gff`。
2. `ref_gap`（肽到最近注释基因的 bp；0 = 重叠，相邻记 1）用来区分"贴边界的延伸"与
   "基因沙漠里的新基因"。**`ref_gap > 2000` 才叫基因沙漠新基因**。

## 5. 参考实测（H99，2026-11-20 定稿）

| 项 | 值 |
|---|---|
| 数据 | 40 raw（CF 4 + TCL 36，37 ℃，TL=Trypsin / LL=LysargiNase）→ **1,837,506 张 MS2**；6 个 RNA-seq 文库（30 ℃ ×3、37 ℃ ×3，合并作表达证据） |
| 参考 | GCF_000149245.1：14 chr + mito，18.9 Mb，7,826 mRNA |
| L0 / LJR | 7,820（与 RefSeq 一致 **7,818/7,820 = 100.0%**）/ 16,492（唯一 14,683） |
| L2 | 一期 150,006 · 二期增量 34,978 · AUGUSTUS 三期 14,267 · 基因间六框 8,000 |
| E1 | **404,069 PSM**，覆盖 7,105/7,820（90.9%） |
| E2 | 一期 1,188,828 / 二期 201,156 / 三期 116,740 / 六框 316 PSM |
| P6 → P6b | 3,951 行 / 3,147 肽 → 候选 3,413 → 合并 **4,016 行**（+65）；Tier 615/903/2,498 |
| **P7 KEEP** | **1,629 行 / 1,006 唯一肽**；DROP 2,387；MISSING **0**；卡值 score ≥ 0.31、cos ≥ 0.90 |
| P8 事件 | 999 肽 / 537 ORF：不同框 379 · N 端 202 · 同框 198 · **基因间 75/42** · C 端 74 · 反义 42 · 内含子 29 · ncRNA 0 |
| 注释 | 537 ORF：strong 20 / moderate 53 / weak 40 / none 424；TM ≥1 **73** |
| 同源 | 16 株 BLAST：物种特异 121 / 复合群 87 / 属 180 / 科 22 / 超出科 48 / 无同源 79（属内 72.3%，属外无同源 38.7%） |
| 交付 | 谱册 1,006 页 + 注释册 537 页 + 案例册 7 页（中英）；结果包各 37 文件；汇报稿 23 页 ×2 |

## 6. 已知坑（真核侧，全部实测）

| 现象 | 真因 | 处理 |
|---|---|---|
| 新肽整层丢失 | 白名单缺 `L2\|` / `LJR\|` | 前缀竖线不是下划线，四个一起查 |
| P6b 报"候选 ORF 0" | `research.POOL_FILES` 是原核池名 | 按真核六池补齐 |
| P6b 任务全崩 `std::bad_alloc` | 注释蛋白 + 候选 ORF 同库 | `--db-mode novel` |
| 四滤形同虚设、MISSING 100 | `_excluders` 用 `L0_`，真核头是 `L0\|` → 排除集空 | 前缀匹配 `L0_` **或** `L0\|` |
| MISSING 复跑后变 KEEP | `--reuse` 把未判定当已判定，走镜像分支 | 只复用 KEEP/DROP |
| 增量卡值被悄悄放宽 | `auto` 在 200 肽子集上重标 Tier1 | 增量必须 `--cards 0.31,0.9` |
| 91% 送检谱静默丢失 | isoPSV Stage A 读 `<run>/unlabeled/` 缓存，`step1_done.txt` 一在就永久跳过 | 换肽集先删 `unlabeled/`（`validate.reset_step1_cache`） |
| 事件全成 ncRNA / intergenic 恒 0 | 基因几何取自 `tx_model.gff` | 见 §4 第 1 条 |
| 内含子肽被标成 N 端延伸 | `classify()` 缺 `intronic` 分支 | 补分支（常量本来就有） |
| AUGUSTUS 模型里的肽全变 ncRNA | 写回 tx_model 缺 CDS 行或 `protein_coding` | 两者都要 |
| 打包后旧图旧报告还在 | 通用 `orfscape pack` 不覆盖 | 用项目 packer（先清空） |

## 7. 工具与脚本

| 脚本 | 作用 |
|---|---|
| `orfscape eukbuild` | 拼接 L0 / LJR / L2（`--preset CREI\|HSAP` 自动找基因组） |
| `orfscape euktx` | 导入/重组转录本，写 `tx_model.gff` |
| `orfscape sframe` | 六框适用性判定与基因间六框建库 |
| `scripts/cneo_l1_augustus.py` | AUGUSTUS gff3 → 转录本 + 增量库 |
| `scripts/cneo_delta_db.py` | 二期/三期增量 ORF 去重建库 |
| `scripts/cneo_homology.py` | 跨物种同源四步（queries / blast / matrix / figs） |
| `scripts/cneo_report.py` / `cneo_cases.py` / `cneo_pdfs.py` / `cneo_pack.py` | 图、案例、册子、结果包 |
