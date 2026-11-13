#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""宿主机侧：盯着 VM 共享盘里的产物标记，文件一落地就取回项目目录。

VM 端（vm_02_l6.sh）每写完一个产物就 `sync` + 原子写 `artifacts/status/<file>.done`；
本脚本轮询这些标记（hgfs 元数据有延迟，用轮询而不是 inotify），校验大小后复制到
`<dest>/<gtf|coverage|merged>/`，并记账到 `<dest>/manifest.tsv`（sha256）。

    # 边跑边取（推荐：先 nohup 起来，再在 VM 里起流程）
    python ORFscape/scripts/vm_watch.py \
        --art  X:/share/_vm_work/l6/artifacts \
        --dest orfscape_projects/<proj>/02_tx --expect 8

    # 只看当前进度
    python ORFscape/scripts/vm_watch.py --art <同上> --dest <同上> --once

退出条件：看到 `ALL_DONE`，或已取得 `--expect` 个 GTF 且 merged 也就位，或超时。
"""
import argparse
import hashlib
import json
import shutil
import sys
import time
from datetime import datetime
from pathlib import Path

KIND_DIR = {"gtf": "gtf", "bedgraph": "coverage", "merged": "merged"}


def _log(msg):
    print("[%s] %s" % (datetime.now().strftime("%H:%M:%S"), msg), flush=True)


def _parse_marker(path):
    out = {}
    try:
        text = path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return None
    for line in text.splitlines():
        if "=" in line:
            k, _, v = line.partition("=")
            out[k.strip()] = v.strip()
    # 兼容 "kind=gtf sample=SRR1" 这种一行多键的写法
    if "kind" in out:
        parts = out["kind"].split()
        out["kind"] = parts[0]
        for extra in parts[1:]:
            if "=" in extra:
                k, _, v = extra.partition("=")
                out.setdefault(k, v)
    return out or None


def _sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _state_path(dest):
    return dest / ".vm_watch_state.json"


def _load_state(dest):
    p = _state_path(dest)
    if p.is_file():
        try:
            return json.loads(p.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            pass
    return {"done": {}}


def _save_state(dest, state):
    _state_path(dest).write_text(json.dumps(state, ensure_ascii=False, indent=1),
                                 encoding="utf-8")


def scan_once(art, dest, expect, state):
    """扫一轮；返回 (新取回数, 是否全部完成)。"""
    st = art / "status"
    if not st.is_dir():
        return 0, False
    got = 0
    all_done = (art / "ALL_DONE").is_file()
    for marker in sorted(st.glob("*.done")):
        key = marker.name
        if key in state["done"]:
            continue
        info = _parse_marker(marker)
        if not info or "file" not in info:
            continue
        name = Path(info["file"]).name
        src = None
        for sub in ("gtf", "cov", "merged"):
            cand = art / sub / name
            if cand.is_file():
                src = cand
                break
        if src is None:
            continue                      # hgfs 还没同步出来，下一轮再说
        try:
            size = src.stat().st_size
        except OSError:
            continue
        want = int(info.get("size") or 0)
        if want and size != want:
            _log("  %s 大小未就绪（%d/%d），稍后再取" % (name, size, want))
            continue
        kind = info.get("kind", "gtf")
        outdir = dest / KIND_DIR.get(kind, "other")
        outdir.mkdir(parents=True, exist_ok=True)
        dst = outdir / name
        shutil.copyfile(src, dst)
        row = "\t".join([datetime.now().isoformat(timespec="seconds"), kind, name,
                         str(size), info.get("lines", ""), _sha256(dst),
                         info.get("sample", ""), info.get("ts", "")])
        man = dest / "manifest.tsv"
        header = "ts\tkind\tfile\tsize\tlines\tsha256\tsample\tvm_ts\n"
        if not man.is_file():
            man.write_text(header, encoding="utf-8")
        with man.open("a", encoding="utf-8") as fh:
            fh.write(row + "\n")
        state["done"][key] = {"file": name, "kind": kind, "size": size}
        _save_state(dest, state)
        got += 1
        _log("取回 %-28s %8d B  (%s)" % (name, size, kind))
    n_gtf = sum(1 for v in state["done"].values() if v.get("kind") == "gtf")
    finished = all_done or (expect > 0 and n_gtf >= expect and any(
        v.get("kind") == "merged" for v in state["done"].values()))
    return got, finished


def main(argv=None):
    ap = argparse.ArgumentParser(description="VM 转录组产物：落地即取")
    ap.add_argument("--art", required=True, help="VM 端 artifacts 目录（共享盘）")
    ap.add_argument("--dest", required=True, help="项目内目标目录（如 02_tx）")
    ap.add_argument("--expect", type=int, default=0, help="期望的样本数（GTF 个数）")
    ap.add_argument("--interval", type=float, default=10.0, help="轮询间隔秒")
    ap.add_argument("--timeout", type=int, default=0, help="超时秒（0=不限）")
    ap.add_argument("--once", action="store_true", help="只扫一轮就退出")
    args = ap.parse_args(argv)

    art = Path(args.art).resolve()
    dest = Path(args.dest).resolve()
    dest.mkdir(parents=True, exist_ok=True)
    if not art.is_dir():
        _log("artifacts 目录还不存在：%s（VM 可能还没开始）" % art)

    state = _load_state(dest)
    t0 = time.time()
    _log("watch %s → %s（已取 %d 个）" % (art, dest, len(state["done"])))
    while True:
        _got, finished = scan_once(art, dest, args.expect, state)
        n = len(state["done"])
        if args.once:
            _log("当前已取 %d 个产物；%s" % (n, "全部完成" if finished else "仍在进行"))
            return 0
        if finished:
            _log("全部到位（%d 个产物），退出" % n)
            _log("下一步：python -m orfscape l6map --workdir <proj> --gff <L6 GFF> --faa <L6 FAA>")
            return 0
        if args.timeout and (time.time() - t0) > args.timeout:
            _log("超时（%ds），已取 %d 个" % (args.timeout, n))
            return 1
        time.sleep(max(args.interval, 1.0))


if __name__ == "__main__":
    sys.exit(main())
