#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""ORFscape 一键跑：P0 → P9，每步之后自动过阶段自检，失败即停。

    python integrations/scripts/run_pipeline.py --workdir orfscape_projects/<proj> --raw <raw 根目录>
    python integrations/scripts/run_pipeline.py ... --from P6 --to P8      # 从中途续跑
    python integrations/scripts/run_pipeline.py ... --skip P6b --dry-run   # 只打印命令

断点续跑：每步有"完成标记"，已存在就跳过（--force 强制重跑）。
日志与 manifest 落在 <workdir>/00_logs/。

包内编排器可用时委派给它（额外写 progress.jsonl 事件流，GUI 可接管）；
否则回落到本脚本自带的阶段循环。
"""
import argparse
import datetime
import json
import os
import shlex
import subprocess
import sys
import time
from pathlib import Path

STEPS = [
    # id,    orfscape 子命令, 自检阶段, 完成标记（glob，相对 workdir）
    ("P0", "scan", "P0", "00_input/raw_info.tsv"),
    ("P1", "genome", "P1", "01_genome/*/genomic.fna"),
    ("D0", "split", "D0", "03_mgf/*.mgf"),
    ("C", "dbbuild", "C", "02_db/*/db_aug.faa"),
    ("E1", "e1", "E1", "04_search_annot/*/*/pFind-Filtered.spectra"),
    ("QC", "qc", "QC", "04a_qc/*"),
    ("E2", "e2", "E2", "04_search_disc/*/*/pFind-Filtered.spectra"),
    ("P6", "novel", "P6", "05_novel/*/novel_peptides.tsv"),
    ("P6b", "research", "P6b", "05b_research/*/novel_after_research.tsv"),
    ("P7", "validate", "P7", "06b_validate_research/novel_validated.tsv"),
]

POST_STEPS = [
    # P8 拆成三个子步，便于单独重跑
    ("P8a", "novgene", "P8", "07_novgene/NG_table.tsv"),
    ("P8b", "novannot", "P8", "08_annotation/*.tsv"),
    ("P8c", "novreport", "P9", "08_report/*.docx"),
]


def find_orfscape_root(start):
    cur = Path(start).resolve()
    for _ in range(6):
        for cand in (cur / "orfscape" / "__init__.py", cur / "ORFscape" / "orfscape" / "__init__.py"):
            if cand.is_file():
                return cand.parent.parent
        if cur.parent == cur:
            break
        cur = cur.parent
    return None


def marker_exists(workdir, pattern):
    return bool(list(Path(workdir).glob(pattern)))


def build_cmd(step_id, sub, args):
    cmd = [sys.executable, "-m", "orfscape", sub, "--workdir", str(Path(args.workdir).resolve())]
    if sub == "scan" and args.raw:
        cmd += ["--raw", str(Path(args.raw).resolve())]
    if sub in ("e1", "e2", "research") and args.thread:
        cmd += ["--thread", str(args.thread)]
    if sub == "e2" and args.dbs:
        cmd += ["--dbs", args.dbs]
    if sub == "research" and args.workers:
        cmd += ["--workers", str(args.workers)]
    if sub == "validate":
        cmd += ["--from-research", "--min-cos", args.min_cos]
        cmd += ["--out", "06b_validate_research"]
        if args.reuse:
            cmd += ["--reuse", args.reuse]
        if args.cards:
            cmd += ["--cards", args.cards]
    if sub in ("novgene", "novannot", "novreport"):
        cmd += ["--from", "06b_validate_research"]
        if args.no_esm and sub == "novannot":
            cmd += ["--no-esm"]
    return cmd


def run_stage_check(workdir, stage):
    script = Path(__file__).resolve().parent / "stage_check.py"
    cmd = [sys.executable, str(script), "--workdir", workdir, "--stage", stage]
    proc = subprocess.run(cmd, capture_output=True, text=True, encoding="utf-8",
                          errors="replace")
    return proc.returncode, (proc.stdout or "") + (proc.stderr or "")


def delegate(argv):
    """优先用包内 orchestrator（写 progress.jsonl 事件流，GUI 可直接接管）。

    深度自检仍用本目录的 stage_check.py：作为 --check-script 传进去。
    """
    for base in Path(__file__).resolve().parents:
        if (base / "orfscape" / "__init__.py").is_file():
            sys.path.insert(0, str(base))
            break
    try:
        from orfscape import orchestrator
    except ImportError:
        return None
    checker = str(Path(__file__).resolve().parent / "stage_check.py")
    args = list(argv)
    if "--check-script" not in args:
        args += ["--check-script", checker]
    return orchestrator.main(args)


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    rc = delegate(argv)
    if rc is not None:
        return rc
    ap = argparse.ArgumentParser(description="ORFscape 一键跑（带阶段自检）")
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--raw", default=None, help="raw 根目录（P0 用）")
    ap.add_argument("--from", dest="start", default="P0", help="起始步（默认 P0）")
    ap.add_argument("--to", dest="end", default="P8c", help="结束步（默认 P8c）")
    ap.add_argument("--skip", default="", help="跳过的步，逗号分隔，如 P6b,QC")
    ap.add_argument("--thread", type=int, default=0, help="pFind 线程数")
    ap.add_argument("--workers", type=int, default=0, help="P6b 并发任务数")
    ap.add_argument("--dbs", default="", help="E2 库后缀（默认工程内设，如 aug,6f）")
    ap.add_argument("--min-cos", default="auto", help="Tier3 相似性卡值（auto 或数值）")
    ap.add_argument("--reuse", default=None, help="P7 增量：复用旧验证目录")
    ap.add_argument("--cards", default=None, help="P7 增量卡值 'score,cos'")
    ap.add_argument("--no-esm", action="store_true", help="novannot 跳过 ESM-2")
    ap.add_argument("--force", action="store_true", help="忽略完成标记，强制重跑")
    ap.add_argument("--dry-run", action="store_true", help="只打印命令")
    args = ap.parse_args(argv)

    workdir = str(Path(args.workdir).resolve())
    if not Path(workdir).is_dir():
        print("项目目录不存在：%s" % workdir)
        return 2
    pkg_root = find_orfscape_root(workdir) or find_orfscape_root(Path(__file__).parent)
    if not pkg_root:
        print("找不到 orfscape 包：请把仓库根加进 PYTHONPATH 或在仓库根运行")
        return 2

    steps = STEPS + POST_STEPS
    ids = [s[0] for s in steps]
    if args.start not in ids or args.end not in ids:
        print("起始/结束步必须是：%s" % "、".join(ids))
        return 2
    lo, hi = ids.index(args.start), ids.index(args.end)
    skip = {s.strip() for s in args.skip.split(",") if s.strip()}

    logs = Path(workdir) / "00_logs"
    logs.mkdir(exist_ok=True)
    stamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    log_path = logs / ("run_%s.log" % stamp)
    manifest = {"workdir": workdir, "raw": args.raw, "started": stamp,
                "pkg_root": str(pkg_root), "steps": []}
    buf = ["ORFscape run %s  workdir=%s" % (stamp, workdir)]

    def flush():
        log_path.write_text("\n".join(buf) + "\n", encoding="utf-8")

    for sid, sub, stage, marker in steps[lo:hi + 1]:
        if sid in skip:
            buf.append("\n-- 跳过 %s" % sid)
            print("跳过 %s" % sid)
            continue
        if not args.force and marker_exists(workdir, marker):
            buf.append("\n-- %s 已完成（%s 存在），跳过" % (sid, marker))
            print("跳过 %s（已完成）" % sid)
            manifest["steps"].append({"step": sid, "status": "skip"})
            continue

        cmd = build_cmd(sid, sub, args)
        line = " ".join(shlex.quote(c) for c in cmd)
        print("\n>>> %s: %s" % (sid, line))
        buf.append("\n>>> %s\n%s" % (sid, line))
        if args.dry_run:
            manifest["steps"].append({"step": sid, "status": "dry-run", "cmd": line})
            continue

        t0 = time.time()
        proc = subprocess.run(cmd, cwd=str(pkg_root), capture_output=True, text=True,
                              encoding="utf-8", errors="replace")
        out = (proc.stdout or "") + (proc.stderr or "")
        buf.append(out[-4000:])
        dt = time.time() - t0
        print("    rc=%d  %.1fs" % (proc.returncode, dt))
        if proc.returncode != 0:
            flush()
            print(out[-2500:])
            print("\n%s 失败（rc=%d），日志：%s" % (sid, proc.returncode, log_path))
            manifest["steps"].append({"step": sid, "status": "fail",
                                      "rc": proc.returncode, "sec": round(dt, 1)})
            (logs / ("manifest_%s.json" % stamp)).write_text(
                json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
            return proc.returncode

        rc, check_out = run_stage_check(workdir, stage)
        buf.append("\n[stage_check %s]\n%s" % (stage, check_out))
        print(check_out.strip())
        if rc != 0:
            flush()
            print("\n%s 产物自检未通过（%s），停在这一步。" % (sid, stage))
            manifest["steps"].append({"step": sid, "status": "check-fail", "sec": round(dt, 1)})
            (logs / ("manifest_%s.json" % stamp)).write_text(
                json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
            return 1
        manifest["steps"].append({"step": sid, "status": "ok", "sec": round(dt, 1)})

    flush()
    (logs / ("manifest_%s.json" % stamp)).write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    print("\n完成。日志：%s" % log_path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
