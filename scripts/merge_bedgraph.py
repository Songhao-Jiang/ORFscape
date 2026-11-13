#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""把多个样本的 bedgraph 合成"平均深度"bedgraph（跨样本按样本数归一）。

    python merge_bedgraph.py --in 02_tx/coverage --out 02_tx/coverage_mean.bedgraph \
        --samples 8 [--gzip-out]

输入：bedtools genomecov -bga 产出的 `chr start end depth`（0-based 半开）。
输出：同格式；`depth = Σ 各样本 depth / --samples`（没覆盖的样本按 0 计），
只写 depth > 0 的区段，保持 0-based 半开。
"""
import argparse
import gzip
import os
import sys
from pathlib import Path


def _open_text(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, encoding="utf-8", errors="replace")


def read_contigs(paths):
    """→ {contig: (length, diff_array)}，两次遍历：先定长度，再累加。"""
    lens = {}
    for p in paths:
        with _open_text(p) as fh:
            for line in fh:
                c = line.split("\t")
                if len(c) < 4:
                    continue
                try:
                    e = int(c[2])
                except ValueError:
                    continue
                if e > lens.get(c[0], 0):
                    lens[c[0]] = e
    diff = {k: [0.0] * (v + 2) for k, v in lens.items()}
    for p in paths:
        with _open_text(p) as fh:
            for line in fh:
                c = line.split("\t")
                if len(c) < 4:
                    continue
                try:
                    s, e, d = int(c[1]), int(c[2]), float(c[3])
                except ValueError:
                    continue
                arr = diff.get(c[0])
                if arr is None or d == 0:
                    continue
                arr[s] += d
                arr[e] -= d
    return lens, diff


def main(argv=None):
    ap = argparse.ArgumentParser(description="多样本 bedgraph → 平均深度")
    ap.add_argument("--in", dest="src", required=True, help="bedgraph 文件或目录")
    ap.add_argument("--out", required=True, help="输出 bedgraph 路径")
    ap.add_argument("--samples", type=int, default=0,
                    help="归一化用的样本数（0=用实际文件数）")
    ap.add_argument("--gzip-out", action="store_true", help="输出 .gz")
    args = ap.parse_args(argv)

    src = Path(args.src).resolve()
    if src.is_dir():
        paths = sorted(p for p in src.iterdir()
                       if p.is_file() and (p.suffix == ".bedgraph" or p.name.endswith(".bedgraph.gz")))
    else:
        paths = [src]
    if not paths:
        print("没有找到 bedgraph：%s" % src)
        return 2
    n = args.samples or len(paths)
    print("合并 %d 个 bedgraph，按 %d 个样本归一" % (len(paths), n))

    lens, diff = read_contigs(paths)
    out_path = Path(args.out).resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if args.gzip_out and not str(out_path).endswith(".gz"):
        out_path = out_path.with_name(out_path.name + ".gz")

    n_seg = 0
    op = gzip.open if str(out_path).endswith(".gz") else open
    with op(out_path, "wt", encoding="utf-8", newline="\n") as fh:
        for contig, arr in sorted(diff.items()):
            cur = 0.0
            start = 0
            val = 0.0
            for i in range(len(arr) - 1):
                cur += arr[i]
                v = round(cur / n, 2)
                if v != val:
                    if val > 0:
                        fh.write("%s\t%d\t%d\t%.2f\n" % (contig, start, i, val))
                        n_seg += 1
                    start, val = i, v
            if val > 0:
                fh.write("%s\t%d\t%d\t%.2f\n" % (contig, start, len(arr) - 1, val))
                n_seg += 1
    size = out_path.stat().st_size
    print("写出 %s（%d 段，%.1f MB）" % (out_path, n_seg, size / 1048576.0))
    return 0


if __name__ == "__main__":
    sys.exit(main())
