#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""预测器基线 CLI 包装（实现见 `orfscape/core/predbench.py`）。

    python scripts/predictor_baseline.py --workdir <proj> [--n 600] [--out <dir>]

等价于 `python -m orfscape predbench --workdir <proj>`；输出写到 `--out`
（缺省 `<proj>/04a_qc`）：predictor_baseline.tsv / predictor_summary.tsv /
predictor_baseline.png。
"""
import argparse
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

from orfscape.core import predbench                 # noqa: E402


def main(argv=None):
    ap = argparse.ArgumentParser(description="预测器基线（注释肽四选一）")
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--n", type=int, default=600)
    ap.add_argument("--out", default=None, help="输出目录（默认 <proj>/04a_qc）")
    args = ap.parse_args(argv)
    res = predbench.run(args.workdir, args.out, n=args.n, log=print)
    print("chosen:", res.get("chosen"))
    return 0 if res.get("chosen") else 1


if __name__ == "__main__":
    raise SystemExit(main())
