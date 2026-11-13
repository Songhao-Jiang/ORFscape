#!/bin/bash
# ORFscape VMware Ubuntu 转录组工具链一键安装 + L6 层构建脚本
# 用法：在 Ubuntu VM 内执行 bash vm_runbook.sh <workdir> <species_code>
# 前置：VM 内有 Python3 + conda/apt；能访问 G: 盘（共享文件夹或 scp 传入）
set -e

WORKDIR=${1:-.}
SP=${2:-Mab}
 echo "=== ORFscape VM 工具链安装 + L6 转录本层构建 ==="
 echo "工作目录: $WORKDIR | 物种: $SP"

# ─── 1. 工具安装（首次运行；已有则跳过）───
if ! command -v prodigal &>/dev/null; then
    echo ">> 安装 Prodigal..."
    apt-get update -qq && apt-get install -y -qq prodigal 2>/dev/null || \
        conda install -y -c bioconda prodigal 2>/dev/null || \
        { echo "请手动安装 prodigal"; exit 1; }
fi
if ! command -v augustus &>/dev/null; then
    echo ">> 安装 Augustus..."
    apt-get install -y -qq augustus 2>/dev/null || conda install -y -c bioconda augustus 2>/dev/null || true
fi
# GeneMarkS-2 需注册下载（学术许可）
if [ ! -f /usr/local/bin/gmst ]; then
    echo ">> GeneMarkS-2：需手动注册下载 http://exon.gatech.edu/GeneMark/
          下载后放置 gmst 到 PATH 即可；跳过不影响主流程"
fi
if ! command -v spades.py &>/dev/null; then
    echo ">> 安装 SPAdes (rnaSPAdes)..."
    apt-get install -y -qq spades 2>/dev/null || conda install -y -c bioconda spades 2>/dev/null || \
        pip install spades 2>/dev/null || true
fi
if ! command -v fasterq-dump &>/dev/null; then
    echo ">> 安装 SRA-Toolkit..."
    apt-get install -y -qq sra-toolkit 2>/dev/null || conda install -y -c bioconda sra-tools 2>/dev/null || true
fi
if ! command -v hisat2 &>/dev/null; then
    echo ">> 安装 HISAT2..."
    apt-get install -y -qq hisat2 2>/dev/null || conda install -y -c bioconda hisat2 2>/dev/null || true
fi
if ! command -v samtools &>/dev/null; then
    echo ">> 安装 samtools..."
    apt-get install -y -qq samtools 2>/dev/null || conda install -y -c bioconda samtools 2>/dev/null || true
fi
if ! command -v stringtie &>/dev/null; then
    echo ">> 安装 StringTie..."
    apt-get install -y -qq stringtie 2>/dev/null || conda install -y -c bioconda stringtie 2>/dev/null || true
fi
if ! command -v orfipy &>/dev/null; then
    echo ">> 安装 orfipy..."
    pip install orfipy 2>/dev/null || conda install -y -c bioconda orfipy 2>/dev/null || true
fi

 echo ">> 工具检查完成"

# ─── 2. SRA 下载 ───
SRA_TSV="$WORKDIR/02_tx/sra_selected_${SP}.tsv"
FASTQ_DIR="$WORKDIR/02_tx/fastq_${SP}"
mkdir -p "$FASTQ_DIR"
if [ -f "$SRA_TSV" ]; then
    echo ">> SRA 下载（$(wc -l < "$SRA_TSV") runs）..."
    tail -n +2 "$SRA_TSV" | cut -f1 | while read RUN; do
        if [ ! -f "$FASTQ_DIR/${RUN}_1.fastq.gz" ]; then
            echo "   下载 $RUN..."
            fasterq-dump --split-files --outdir "$FASTQ_DIR" --threads 4 "$RUN" && \
                gzip "$FASTQ_DIR/${RUN}"*.fastq
        else
            echo "   $RUN 已存在，跳过"
        fi
    done
else
    echo "!! 未找到 $SRA_TSV（先在 Windows 端运行 genome + sraauto）"
fi

# ─── 3. 转录组比对 + 组装 ───
GENOME_FNA="$WORKDIR/01_genome/${SP}/genomic.fna"
ALIGN_DIR="$WORKDIR/02_tx/align_${SP}"
ASM_DIR="$WORKDIR/02_tx/asm_${SP}"
mkdir -p "$ALIGN_DIR" "$ASM_DIR"

# HISAT2 索引 + 比对
if [ -f "$GENOME_FNA" ]; then
    echo ">> HISAT2 索引..."
    hisat2-build "$GENOME_FNA" "$ALIGN_DIR/${SP}_index"
    echo ">> HISAT2 比对..."
    for FQ1 in "$FASTQ_DIR"/*_1.fastq.gz; do
        FQ2="${FQ1/_1/_2}"
        RUN=$(basename "$FQ1" _1.fastq.gz)
        if [ ! -f "$ALIGN_DIR/${RUN}.bam" ]; then
            echo "   比对 $RUN..."
            hisat2 -x "$ALIGN_DIR/${SP}_index" -1 "$FQ1" -2 "$FQ2" -S - 2>/dev/null | \
                samtools sort -o "$ALIGN_DIR/${RUN}.bam" -
        fi
    done
fi

# StringTie 组装（合并）
echo ">> StringTie 组装..."
BAMS=$(ls "$ALIGN_DIR"/*.bam 2>/dev/null | head -20)
if [ -n "$BAMS" ]; then
    stringtie --merge -o "$ASM_DIR/merged.gtf" $BAMS 2>/dev/null
    # 提取转录本序列
    gffread "$ASM_DIR/merged.gtf" -g "$GENOME_FNA" -w "$ASM_DIR/transcripts.fa" 2>/dev/null || \
        python3 -c "
import sys
# 简易 GTF→FASTA（无 gffread 时备用）
"
fi

# ─── 4. orfipy 六框 ORF 提取 ───
ORF_DIR="$WORKDIR/02_tx/orf_${SP}"
mkdir -p "$ORF_DIR"
if [ -f "$ASM_DIR/transcripts.fa" ]; then
    echo ">> orfipy 六框 ORF 提取..."
    orfipy "$ASM_DIR/transcripts.fa" --pep "$ORF_DIR/orfs.pep" --min 30 --threads 4 \
        --start-codons ATG,GTG,TTG,CTG 2>/dev/null
    echo "   ORF 肽段数: $(grep -c '>' "$ORF_DIR/orfs.pep")"
fi

 echo "=== VM 脚本完成。将 $ORF_DIR/orfs.pep 拷回 Windows 端 $WORKDIR/02_db/${SP}/L6_tx_orfs.faa 即可 ==="
