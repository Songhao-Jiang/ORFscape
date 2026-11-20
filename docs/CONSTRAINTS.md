# Constraints — 限定文本

**English summary.** These are the rules that make numbers comparable and runs
reproducible. They are not style preferences: changing any of them changes the
result, and a result computed under a different rule must never be reported next to
one computed under this document. Each rule states *why*, so a future maintainer can
decide whether to re-baseline instead of silently drifting.

> 用法：写方法学/报告/答辩时，凡引用 KEEP、事件、ORF、注释数，都默认受本文约束；
> 改动任一条 → 视为新口径，必须重跑并在记录里写明"与旧版不可比"。

## 1. 数字可比性（最容易被破坏的三条）

| # | 规则 | 为什么 |
|---|---|---|
| 1.1 | **卡值（score / cos）每轮按 Tier1 重新标定，跨轮不可比**；报告 KEEP 时必须同时给卡值 | 不同数据的 Tier1 分布不同，0.82 与 0.86 的 KEEP 不能直接比大小 |
| 1.2 | 增量模式**禁止 `--cards auto`**（子集只有几百肽，Tier1 锚定会把门控放宽） | 实测 0.31/0.90 被 auto 拉到 0.23/0.70（n=11），整表门控失效 |
| 1.3 | **`--reuse` 只复用 KEEP / DROP**，MISSING 一律重验 | 把"未判定"当"已判定"会经镜像分支提升成 KEEP（实测 +99 肽） |

## 2. 建库与层

| # | 规则 | 为什么 |
|---|---|---|
| 2.1 | 层前缀必须完整覆盖：原核 `L1_` / `L2_` / `L3\|` / `L6\|`；真核加 `LJR\|` / `L2\|` / `L1_SF\|` | 漏一个 → 整层静默丢失，日志不报错 |
| 2.2 | 序列一律**自己从基因组 + 坐标重拼**，不信厂家 FASTA | FASTA 与坐标必须严格自洽；参考校验是唯一硬证据 |
| 2.3 | 翻译读 GFF 的 `transl_table`（含线粒体/霉菌线粒体表） | 写死表 1 会把线粒体与部分真菌基因翻错 |
| 2.4 | 真核 AUGUSTUS 模型写回 `tx_model.gff` 必须带 **CDS 行 + `protein_coding`** | 否则整层被判成非编码转录本，肽全记成 `ncRNA_translation` |
| 2.5 | 六框只在**通过 R 门**后建；整库与基因间分建，先看基因间占比 | R = 6/d；稀疏基因组做六框只会放大库、压低真肽分数 |
| 2.6 | 诱饵与目标等量（`REV_` 前缀），由 pFind 内部估 FDR | 自己算 FDR 会与搜索引擎口径冲突 |
| 2.7 | L6 必须在 `dbbuild` **之前**完成；bedgraph 必须在删 BAM 之前导出 | 事后补不回来；merged.gtf 没有丰度 |

## 3. 搜库

| # | 规则 | 为什么 |
|---|---|---|
| 3.1 | 全项目同一套 pFind 参数（半特异 `digest=1`、20 ppm、Carbamidomethyl(C)、Oxidation(M)/Acetyl(N-term)、漏切 2、I/L 等价、PSM FDR 1%） | 参数一变，E1/E2/P6b 之间不可比 |
| 3.2 | E2 按库分池并行，**互不混池**（原核 aug / 6f；真核 6f / 6g / 8a / 9i） | 报告里要能分开说"转录证据"与"从头预测" |
| 3.3 | 结果目录必须带**库指纹**（`db.sha256`），库变则旧结果作废 | 否则换库后静默复用旧命中 |
| 3.4 | P6b 回搜库**只放候选 ORF**（`--db-mode novel`） | 注释蛋白同库 → pFind `Group.exe` 组合爆炸（`std::bad_alloc`），40/40 任务崩 |

## 4. 新肽与回搜

| # | 规则 | 为什么 |
|---|---|---|
| 4.1 | 四滤全过才算候选：① 不在 L0 双酶消化集 ② 不在 E1 实测肽集 ③ 不是 L0 蛋白 ≥7 aa 连续子串 ④ 等重变体同样要过 ①② | 任一环节空集，等于四滤失效（真核实测漏 100 条已注释肽） |
| 4.2 | I/L 归一是**所有**比较的前提（库、消化集、变体） | 质谱分不开 I/L；不归一会出现"新肽"其实是已知肽 |
| 4.3 | 合并集行数必须 ≥ P6；否则视为回搜无效 | 回搜只会加证据，不会减肽 |
| 4.4 | 回搜失败率 ≥ 50% 直接抛错，不写"0 命中" | 0 命中是故障信号，不是结果 |

## 5. 判读（P7）

| # | 规则 | 为什么 |
|---|---|---|
| 5.1 | 每条送检肽都必须有结论；`MISSING` 只能表示"未判定"，且必须带原因 | 未判读者不得计入保留 |
| 5.2 | 换肽集必须清空 isoPSV 的 `<run_dir>/unlabeled/`（`step1_done.txt` 同删） | Stage A 读缓存，不清就静默丢谱（实测 91%） |
| 5.3 | 已注释库命中者按 `proteogenomics_mode=1` 跳过错误校正，记 `annotated_skipped` | 这类肽随后被四滤剔除，不进结果表 |
| 5.4 | Tier1 = 同肽跨酶；Tier2 = 同 locus 双酶；Tier3 = 单酶。Tier3 必须"三门全过" | 单酶证据最弱，不能与跨酶同权 |
| 5.5 | Tier 变动（T3→T2 等）的肽**不能用旧结论**，必须按新阈值重验 | 卡值按 Tier 标定，跨档后旧阈值不成立 |
| 5.6 | isoPSV 版本必须写进汇总表头（由安装路径自动取） | 曾经硬编码 v2.0.175，实际跑 2.0.184 |

## 6. 事件与计数

| # | 规则 | 为什么 |
|---|---|---|
| 6.1 | 行 = 肽 × 来源 ORF；对外报**唯一肽**（I/L 归一后去重）。两个数不许混用 | 一行不等于一条肽 |
| 6.2 | 事件按**肽级归一**：一条肽只记一个主事件 | 否则同一肽在多 ORF 上重复计数 |
| 6.3 | 原核事件互斥优先级：**N 端校正 > 内部不同框 > 反义 > 基因间** | 早发现的事件信息量更高 |
| 6.4 | 事件只比**同链** CDS；重叠多个 CDS 取 `rel % 3 == 0` 的同框参照 | 否则 N 端事件虚高 |
| 6.5 | ORF 去重：同 contig / 同链 / **同读框** / 重叠 > 80% 视为同一位点 | 别名 ORF 会把 ORF 数翻倍 |
| 6.6 | 真核基因几何**只用参考注释** `01_genome/<sp>/genomic.gff`；坐标映射用 `tx_model.gff` | 用 tx_model 会把 span 撑进基因间区（intergenic 恒 0） |
| 6.7 | `ref_gap > 2000` 才称"基因沙漠新基因"；`≤ 500` 多半是注释边界不准 | 贴边界事件不是新基因 |
| 6.8 | 同源特异性阈值：蛋白 id ≥ 50 且 cov ≥ 70，或 核酸 id ≥ 70 且 cov ≥ 70；cov 用查询跨度 `qend−qstart+1` | 用 `length/qlen` 会把 gapped HSP 算小，计数漂移 |

## 7. 工程约束

| # | 规则 | 为什么 |
|---|---|---|
| 7.1 | 所有写盘先校验归属：`abspath` + `commonpath == base`（模块内 `_open_w(base, path, mode)`） | 防止 `--out` 相对路径把产物写进别的工程 |
| 7.2 | `--out` 一律绝对路径 | 相对路径会落到 ORFscape 仓库根 |
| 7.3 | 外部输入不做字符串拼接 SQL；URL 只允许 http/https 且校验 host（拒绝 localhost/环回/私网/保留地址、禁重定向） | 安全基线 |
| 7.4 | 解析不可信 XML 前拒绝 `DOCTYPE`/`ENTITY`，不开外部实体 | 安全基线 |
| 7.5 | 凭据只从环境变量或密钥服务读；源码、示例、测试不写可用凭据 | 安全基线 |
| 7.6 | 不新增依赖（GUI 只用 tkinter + 标准库） | 保证离线可跑 |
| 7.7 | 结果包必须**先清空再拷**；不复用旧目录 | 通用 `orfscape pack` 不覆盖，会留下上一轮的图与报告 |

## 8. 反模式速查

| 写法 | 后果 | 正确写法 |
|---|---|---|
| `"L" in "TL"` | LL 被当成 TL，酶判定全错 | `== "TL"` / `== "LL"` |
| `h.startswith("L0_")`（真核） | 四滤排除集为空 | `startswith("L0_") or startswith("L0\|")` |
| `--min-cos auto` 跑增量 | 门控被放宽 | `--cards 0.31,0.9` |
| `validate` 前不删 `unlabeled/` | 新谱不进判读 | `reset_step1_cache` |
| `pack` 直接跑到旧目录 | 包内数字与图不一致 | 先 `rm -rf` 目标 |
| 事件基因集用 `tx_model.gff` | intergenic 恒 0 | 参考 GFF |
| `makeblastdb` 建到网络盘 | 无法创建目录 | 本地盘 |
| 报告 KEEP 498 / 934 / 1,006 而不说卡值 | 读者以为是同口径 | 数字 + 卡值 + 版本一起给 |

## 9. 改口径的正确姿势

1. 在实验记录里新开一节，写"改了什么、为什么、由谁定"；
2. 全量重跑受影响的阶段（不是只跑增量）；
3. 新表标题写"新口径"，旧数字保留但标注"作废，不可并列"；
4. 同步更新方法学与交接文件里的预期数字与门；
5. 若改的是代码而非参数，先跑 `tests/` 全量，再跑一次幂等检查（连续两次结果一致）。
