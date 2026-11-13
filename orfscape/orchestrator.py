# -*- coding: utf-8 -*-
"""阶段编排：P0 → P9，带流式日志、结构化事件与产物自检。

与 ``.zcode/skills/orfscape-proteogenomics/scripts/run_pipeline.py`` 同语义，
但自包含（不依赖 skill 目录），并额外做三件事：

1. 子进程输出逐行转发，同时由子进程自己写 ``ORFSCAPE_EVENTS`` 事件流；
2. 每个阶段边界补 ``stage.start / stage.end / stage.check`` 事件，GUI 据此
   画看板与 ETA；
3. 结束（含失败/停止）时把实测阶段耗时写进 ``00_logs/calibration.json``，
   下一次预估直接用真值。

停止：GUI 写 ``00_logs/stop.flag``，本模块在行边界检测到就终止当前子进程并
优雅收尾。

路径安全：所有项目内路径都由"项目根 + 字面量子目录 + 单段文件名"拼成，
文件名先过白名单正则，再校验 dirname 必须等于基目录。
"""
import argparse
import json
import os
import re
import subprocess
import sys
import time
from pathlib import Path

from . import checks, estimator, progress
from .cli import ensure_project

# id, 子命令, 自检阶段, 完成标记(glob), 权重
STEPS = [
    ("P0", "scan", "P0", "00_input/raw_info.tsv", 1.0),
    ("P1", "genome", "P1", "01_genome/*/genomic.fna", 3.0),
    ("D0", "split", "D0", "03_mgf/*.mgf", 10.0),
    ("C", "dbbuild", "C", "02_db/*/db_aug.faa", 2.0),
    ("E1", "e1", "E1", "04_search_annot/*/*/pFind-Filtered.spectra", 9.0),
    ("QC", "qc", "QC", "04a_qc/*", 0.4),
    ("E2", "e2", "E2", "04_search_disc/*/*/pFind-Filtered.spectra", 30.0),
    ("P6", "novel", "P6", "05_novel/*/novel_peptides.tsv", 1.5),
    ("P6b", "research", "P6b", "05b_research/*/novel_after_research.tsv", 14.0),
    ("P7", "validate", "P7", "06b_validate_research/novel_validated.tsv", 12.0),
    ("P8a", "novgene", "P8", "07_novgene/NG_table.tsv", 0.6),
    ("P8b", "novannot", "P8", "08_annotation/*.tsv", 4.0),
    ("P8c", "novreport", "P9", "08_report/*.docx", 2.0),
]

TITLES = {
    "P0": "Data confirm", "P1": "Genome fetch", "D0": "Spectra / CV split",
    "C": "Layer database", "E1": "Annotated search", "QC": "QC",
    "E2": "Discovery search", "P6": "Novel peptides", "P6b": "Candidate re-search",
    "P7": "isoPSV arbitration", "P8a": "Event classification",
    "P8b": "Functional annotation", "P8c": "Report",
}

CHECK_TITLES = {
    "P0": "raw list", "P1": "assembly", "D0": "MGF", "C": "layer FASTA",
    "E1": "annotated PSM", "QC": "qc outputs", "E2": "discovery PSM",
    "P6": "novel peptides", "P6b": "re-search peptides", "P7": "validated",
    "P8": "events", "P9": "report",
}

PARENT = os.pardir
_PART = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._\-]*$")


def _now():
    return time.strftime("%H:%M:%S")


# ---------------------------------------------------------------- 路径安全
def project_root(workdir):
    """项目根：拒绝父目录成分 + 规范化 + 必须已存在。"""
    raw = str(workdir).replace("\\", "/")
    if any(part == PARENT for part in raw.split("/")):
        raise ValueError("project path contains a parent component")
    root = os.path.realpath(os.path.abspath(str(workdir)))
    if not os.path.isdir(root):
        raise ValueError("project directory not found")
    return root


def checked_join(base, name):
    """base 下拼一个**单段**文件名：白名单 + dirname 必须等于 base。"""
    name = str(name)
    if not _PART.match(name) or name in (".", PARENT):
        raise ValueError("unsafe file name: %s" % name)
    base_r = os.path.realpath(base)
    out = os.path.realpath(os.path.join(base_r, name))
    if os.path.dirname(out) != base_r:
        raise ValueError("file escapes its directory: %s" % name)
    return out


def rm_checked(path, root):
    """只删 root 内、真实存在的普通文件。"""
    root_r = os.path.realpath(root)
    p = os.path.realpath(path)
    if os.path.commonpath([p, root_r]) != root_r:
        raise ValueError("refuse to remove outside project: %s" % p)
    if os.path.isfile(p):
        os.remove(p)
        return True
    return False


def marker_exists(workdir, pattern):
    """项目根下的相对 glob（禁父目录段，命中项必须在项目内）。"""
    rel = str(pattern).replace("\\", "/")
    segs = [s for s in rel.split("/") if s not in ("", ".")]
    if not segs or any(s == PARENT for s in segs):
        raise ValueError("unsafe marker pattern")
    base = Path(os.path.realpath(workdir))
    for hit in base.glob(rel):
        if hit == base or base in hit.parents:
            return True
    return False


def input_path(path):
    """用户给的 raw 路径：拒绝父目录成分后规范化（允许在项目外）。"""
    raw = str(path).replace("\\", "/")
    if any(part == PARENT for part in raw.split("/")):
        raise ValueError("input path contains a parent component")
    return os.path.realpath(os.path.abspath(str(path)))


# ---------------------------------------------------------------- 阶段表
def find_pkg_root(start):
    """向上最多 7 层找 orfscape 包。"""
    cur = Path(os.path.realpath(os.path.abspath(str(start))))
    for _ in range(7):
        if (cur / "orfscape" / "__init__.py").is_file():
            return str(cur)
        if (cur / "ORFscape" / "orfscape" / "__init__.py").is_file():
            return str(cur / "ORFscape")
        if cur.parent == cur:
            break
        cur = cur.parent
    return None


def build_cmd(step_id, sub, args, workdir):
    cmd = [sys.executable, "-m", "orfscape", sub, "--workdir", workdir]
    if sub == "scan":
        if args.raw:
            cmd += ["--raw", args.raw]
        if getattr(args, "species_name", None):
            cmd += ["--species-name", args.species_name]
        if getattr(args, "species_code", None):
            cmd += ["--species-code", args.species_code]
    if sub in ("e1", "e2") and args.thread:
        cmd += ["--thread", str(args.thread)]
    if sub == "e2" and args.dbs:
        cmd += ["--dbs", args.dbs]
    if sub == "research" and args.workers:
        cmd += ["--workers", str(args.workers)]
    if sub == "validate":
        cmd += ["--from-research", "--min-cos", args.min_cos,
                "--out", "06b_validate_research"]
        if args.reuse:
            cmd += ["--reuse", args.reuse]
        if args.cards:
            cmd += ["--cards", args.cards]
    if sub in ("novgene", "novannot", "novreport"):
        cmd += ["--from", "06b_validate_research"]
        if args.no_esm and sub == "novannot":
            cmd += ["--no-esm"]
    return cmd


def run_check(workdir, stage, script=None):
    """阶段自检：优先用外部脚本（skill 的深度检查），否则用内置 checks。"""
    if not script:
        return checks.check_stage(workdir, stage)
    cmd = [sys.executable, script, "--workdir", workdir, "--stage", stage]
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True,
                              encoding="utf-8", errors="replace")
    except OSError as exc:
        return {"status": "warn", "fails": 0, "warns": 1,
                "lines": ["  WARN checker failed to start: %s" % exc],
                "summary": "checker unavailable"}
    out = ((proc.stdout or "") + (proc.stderr or "")).splitlines()
    fails = warns = 0
    for line in out:
        if "FAIL" in line and "/" in line and ("WARN" in line or "合计" in line):
            parts = line.replace("合计：", "").replace("合计:", "").split("/")
            try:
                fails = int("".join(ch for ch in parts[0] if ch.isdigit()) or 0)
                warns = int("".join(ch for ch in parts[1] if ch.isdigit()) or 0)
                break
            except (IndexError, ValueError):
                continue
    summary = out[-1] if out else ""
    status = "fail" if (fails or proc.returncode) else ("warn" if warns else "ok")
    return {"stage": stage, "status": status, "fails": fails, "warns": warns,
            "lines": out, "summary": summary}


def project_paths(workdir):
    """一次算好本轮要用的全部路径（都在项目内，字面量子目录）。"""
    root = project_root(workdir)
    logs = os.path.join(root, "00_logs")
    os.makedirs(logs, exist_ok=True)
    stamp = time.strftime("%Y%m%d_%H%M%S")
    return {
        "root": root,
        "logs": logs,
        "events": checked_join(logs, "progress.jsonl"),
        "stop": checked_join(logs, "stop.flag"),
        "log": checked_join(logs, "run_%s.log" % stamp),
        "manifest": checked_join(logs, "manifest_%s.json" % stamp),
        "stamp": stamp,
    }


# ---------------------------------------------------------------- 主流程
def run(args):
    """执行整条流程；返回进程退出码。"""
    try:
        pth = project_paths(args.workdir)
    except (ValueError, OSError) as exc:
        print("bad project path: %s" % exc)
        return 2
    workdir = pth["root"]

    pkg_root = find_pkg_root(workdir) or find_pkg_root(os.path.dirname(__file__))
    if not pkg_root:
        print("cannot locate the orfscape package; run `pip install -e ORFscape`")
        return 2

    if os.path.isfile(pth["stop"]):
        rm_checked(pth["stop"], workdir)
    # 新一轮：重建事件文件，避免 GUI 读到上一轮的行
    with open(pth["events"], "w", encoding="utf-8"):
        pass
    os.environ[progress.ENV_PATH] = pth["events"]
    if args.lang:
        os.environ["ORFSCAPE_LANG"] = args.lang

    log_fh = open(pth["log"], "w", encoding="utf-8", newline="\n")

    def say(text, event=True, stage=None):
        line = "%s %s" % (_now(), text)
        print(line, flush=True)
        log_fh.write(line + "\n")
        log_fh.flush()
        if event:
            progress.emit("log", line=line, stage=stage)

    ids = [s[0] for s in STEPS]
    if args.start not in ids or args.end not in ids:
        print("--from/--to must be one of: %s" % ", ".join(ids))
        return 2
    lo, hi = ids.index(args.start), ids.index(args.end)
    skip = {s.strip() for s in (args.skip or "").split(",") if s.strip()}

    est = estimator.estimate(workdir, [args.raw] if args.raw else None,
                             threads=args.thread, workers=args.workers,
                             dbs=args.dbs, species_code=args.species_code)

    selected = list(STEPS[lo:hi + 1])
    progress.emit("run.start", workdir=workdir, raw=args.raw,
                  species=args.species_name or args.species_code,
                  desc=args.desc, lang=args.lang or "en",
                  stages=[{"id": sid, "title_key": "stage." + sid,
                           "title": TITLES.get(sid, sid), "weight": w,
                           "est_sec": (est["est_sec"].get(sid) or {}).get("sec")}
                          for sid, _sub, _chk, _mk, w in selected],
                  estimate={"total_sec": est["total_sec"],
                            "disk_bytes": est["disk"]["total"],
                            "ram_gb": est["ram"]["need_gb"],
                            "ms2": est["inputs"]["ms2"],
                            "n_raw": est["inputs"]["n_raw"]})
    say("ORFscape run %s  project=%s" % (pth["stamp"], workdir))

    manifest = {"workdir": workdir, "raw": args.raw, "started": pth["stamp"],
                "pkg_root": pkg_root, "steps": []}
    started_all = time.time()
    rc_all = 0
    status_all = "ok"

    try:
        for sid, sub, chk, marker, weight in selected:
            if sid in skip:
                say("-- skip %s" % sid, stage=sid)
                progress.emit("stage.start", stage=sid, skipped=True)
                progress.emit("stage.end", stage=sid, status="skip", sec=0)
                manifest["steps"].append({"step": sid, "status": "skip"})
                continue
            if not args.force and marker_exists(workdir, marker):
                say("-- %s already done, skip" % sid, stage=sid)
                progress.emit("stage.start", stage=sid, skipped=True)
                progress.emit("stage.end", stage=sid, status="skip", sec=0)
                manifest["steps"].append({"step": sid, "status": "skip"})
                continue

            cmd = build_cmd(sid, sub, args, workdir)
            say(">>> %s: %s" % (sid, " ".join(cmd)), stage=sid)
            if args.dry_run:
                manifest["steps"].append({"step": sid, "status": "dry-run",
                                          "cmd": " ".join(cmd)})
                continue

            progress.emit("stage.start", stage=sid, cmd=" ".join(cmd),
                          weight=weight,
                          est_sec=(est["est_sec"].get(sid) or {}).get("sec"),
                          check=CHECK_TITLES.get(chk, chk))
            t0 = time.time()
            p_rc = _run_child(cmd, pkg_root, pth["stop"], log_fh)
            dt = time.time() - t0

            if p_rc is None:                     # 被 stop.flag 打断
                status_all = "stopped"
                progress.emit("stage.end", stage=sid, status="stopped",
                              sec=round(dt, 1))
                say("-- stopped by user at %s" % sid)
                manifest["steps"].append({"step": sid, "status": "stopped",
                                          "sec": round(dt, 1)})
                rc_all = 130
                break

            say("    rc=%d  %.1fs" % (p_rc, dt), stage=sid)
            if p_rc != 0:
                status_all = "fail"
                progress.emit("stage.end", stage=sid, status="fail",
                              sec=round(dt, 1), rc=p_rc)
                manifest["steps"].append({"step": sid, "status": "fail",
                                          "rc": p_rc, "sec": round(dt, 1)})
                rc_all = p_rc or 1
                break

            res = run_check(workdir, chk, getattr(args, "check_script", None))
            say("[stage_check %s] %s" % (chk, res.get("summary") or ""),
                stage=sid)
            progress.emit("stage.check", stage=sid, status=res.get("status"),
                          fails=res.get("fails"), warns=res.get("warns"),
                          lines=res.get("lines") or [])
            if res.get("status") == "fail":
                status_all = "check-fail"
                progress.emit("stage.end", stage=sid, status="check-fail",
                              sec=round(dt, 1))
                manifest["steps"].append({"step": sid, "status": "check-fail",
                                          "sec": round(dt, 1)})
                rc_all = 1
                break

            progress.emit("stage.end", stage=sid, status="ok", sec=round(dt, 1))
            manifest["steps"].append({"step": sid, "status": "ok",
                                      "sec": round(dt, 1)})
    except Exception as exc:                 # 编排器自身异常也要收尾
        status_all = "fail"
        rc_all = rc_all or 1
        say("orchestrator error: %r" % exc, event=False)
        progress.emit("error", text="%r" % exc)
    finally:
        total = time.time() - started_all
        progress.emit("run.end", status=status_all, sec=round(total, 1),
                      rc=rc_all)
        manifest["ended"] = time.strftime("%Y%m%d_%H%M%S")
        manifest["sec"] = round(total, 1)
        manifest["status"] = status_all
        try:
            with open(pth["manifest"], "w", encoding="utf-8") as fh:
                json.dump(manifest, fh, ensure_ascii=False, indent=2)
        except OSError:
            pass
        try:
            from .runlog import RunState, read_events
            st = RunState()
            for rec in read_events(pth["events"])[0]:
                st.feed(rec)
            if estimator.save_calibration(workdir, st):
                say("calibration saved", event=False)
        except Exception:
            pass
        say("done. log: %s" % pth["log"], event=False)
        log_fh.close()
    return rc_all


def _run_child(cmd, cwd, stop_flag, log_fh):
    """跑一个子命令，逐行转发；收到 stop.flag 返回 None。"""
    env = dict(os.environ)
    env["PYTHONUNBUFFERED"] = "1"
    env["PYTHONIOENCODING"] = "utf-8"
    try:
        proc = subprocess.Popen(cmd, cwd=cwd, env=env, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, text=True,
                                encoding="utf-8", errors="replace", bufsize=1)
    except OSError as exc:
        print("spawn failed: %s" % exc, flush=True)
        return 127
    stopped = False
    try:
        assert proc.stdout is not None
        for line in proc.stdout:
            line = line.rstrip("\n")
            print(line, flush=True)
            log_fh.write(line + "\n")
            log_fh.flush()
            if os.path.isfile(stop_flag):
                stopped = True
                break
    finally:
        if stopped:
            _kill_tree(proc)
            try:
                proc.wait(timeout=20)
            except subprocess.TimeoutExpired:
                pass
            try:
                rm_checked(stop_flag, os.path.dirname(stop_flag))
            except (OSError, ValueError):
                pass
            return None
        rc = proc.wait()
    return rc


def _kill_tree(proc):
    """终止子进程树（python 再拉起 pFind 子进程）。"""
    try:
        if os.name == "nt":
            subprocess.run(["taskkill", "/F", "/T", "/PID", str(proc.pid)],
                           capture_output=True)
        else:
            proc.terminate()
    except OSError:
        try:
            proc.kill()
        except OSError:
            pass


def main(argv=None):
    ap = argparse.ArgumentParser(description="ORFscape staged runner (events + checks)")
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--raw", default=None)
    ap.add_argument("--species-name", default=None,
                    help="species Latin name; overrides filename inference")
    ap.add_argument("--species-code", default=None,
                    help="species code, e.g. TP (default: inferred)")
    ap.add_argument("--desc", default=None, help="free-text description")
    ap.add_argument("--from", dest="start", default="P0")
    ap.add_argument("--to", dest="end", default="P8c")
    ap.add_argument("--skip", default="")
    ap.add_argument("--thread", type=int, default=0)
    ap.add_argument("--workers", type=int, default=0)
    ap.add_argument("--dbs", default="aug,6f")
    ap.add_argument("--min-cos", default="auto")
    ap.add_argument("--reuse", default=None)
    ap.add_argument("--cards", default=None)
    ap.add_argument("--no-esm", action="store_true")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--check-script", default=None,
                    help="external stage checker (e.g. the skill's stage_check.py)")
    ap.add_argument("--lang", default=None, choices=["en", "zh"])
    args = ap.parse_args(argv)
    if not args.thread:
        args.thread = max(4, (os.cpu_count() or 8) // 2)
    if not args.workers:
        args.workers = 4
    if args.raw:
        args.raw = input_path(args.raw)
    if not args.dry_run:
        ensure_project(args.workdir)
    return run(args)


if __name__ == "__main__":
    sys.exit(main())
