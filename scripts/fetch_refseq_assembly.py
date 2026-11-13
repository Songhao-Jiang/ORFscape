# -*- coding: utf-8 -*-
"""从 NCBI RefSeq 下载一个装配的三件套（genomic.fna / genomic.gff / protein.faa）。

真核演示用：原核走 `orfscape genome`（按物种名查 assembly_summary），这里直接
按装配号钉扎下载，避免物种名解析的歧义。

    python ORFscape/scripts/fetch_refseq_assembly.py \
        --stem GCF_000146045.2_R64 --out orfscape_projects/euk_pilot/01_genome/SCER

安全：只允许 https + ftp.ncbi.nlm.nih.gov（复用 genome._guarded_urlopen 的校验），
不跟随重定向；装配号与 stem 都过白名单正则，路径只用单段文件名。
"""
import argparse
import gzip
import os
import shutil
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

from orfscape.core import genome as G                    # noqa: E402

_RE_STEM = re.compile(r"^GC[AF]_\d+\.\d+_[A-Za-z0-9_.\-]+$")
_SUFFIX = {"genomic.fna": "_genomic.fna.gz",
           "genomic.gff": "_genomic.gff.gz",
           "protein.faa": "_protein.faa.gz"}


def url_for(stem, kind):
    acc = "_".join(stem.split("_")[:2])          # GCF_000146045.2
    parts = acc.split("_")[-1].split(".")[0]     # 000146045
    if len(parts) < 9:
        return None
    d1, d2, d3 = parts[0:3], parts[3:6], parts[6:9]
    return ("https://ftp.ncbi.nlm.nih.gov/genomes/all/%s/%s/%s/%s/%s/%s%s"
            % (acc.split("_")[0], d1, d2, d3, stem, stem, _SUFFIX[kind]))
    # 例：.../all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stem", required=True, help="完整装配名，如 GCF_000146045.2_R64")
    ap.add_argument("--out", required=True)
    ap.add_argument("--files", default="genomic.fna,genomic.gff,protein.faa")
    args = ap.parse_args()

    if not _RE_STEM.match(args.stem):
        print("非法装配名：%s" % args.stem)
        return 2
    out = Path(os.path.realpath(args.out))
    out.mkdir(parents=True, exist_ok=True)

    for kind in [x.strip() for x in args.files.split(",") if x.strip()]:
        if kind not in _SUFFIX:
            print("跳过未知文件：%s" % kind)
            continue
        url = url_for(args.stem, kind)
        if not url:
            print("装配号格式不对：%s" % args.stem)
            return 2
        target = out / kind
        print("下载 %s" % url, flush=True)
        # 流式解压落盘：人基因组 gz 250 MB / 解压 3.1 GB，整块读进内存没必要
        with G._guarded_urlopen(url) as resp, gzip.GzipFile(fileobj=resp) as gz:
            with target.open("wb") as fh:
                shutil.copyfileobj(gz, fh, 1 << 20)
        print("  → %s（%.1f MB）" % (target, target.stat().st_size / 1e6),
              flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
