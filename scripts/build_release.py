#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Build the ORFscape release asset: a clean tree -> dist/ORFscape_v<ver>.zip.

    python scripts/build_release.py              # build + verify
    python scripts/build_release.py --check      # scan only, write nothing

What it does
  1. walks the repository, skipping caches, local config, datasets and internal notes;
  2. scans every text file for personal paths / credentials and **fails** the build
     if one is found (a release must not contain them);
  3. writes dist/ORFscape_v<ver>.zip with a single top-level `ORFscape/` directory;
  4. prints exactly one SHA-256, computed after the archive is closed, and records it
     in dist/SHA256SUMS.txt (no arithmetic the reader has to redo).
"""
import argparse
import hashlib
import os
import re
import sys
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DIST = ROOT / "dist"

# 目录整棵跳过
SKIP_DIRS = {
    "__pycache__", ".git", ".mimosa", ".idea", ".vscode", "dist", "build",
    "orfscape_projects", "raw", "00_logs", ".venv", "venv", "node_modules",
    "orfipy_transcribed_regions.fa_out",
}
# 文件跳过（相对仓库根，glob 也可）
SKIP_GLOBS = [
    "*.pyc", "*.pyo", "*.log", "*.jsonl", "*.zip", "*.tmp",
    "orfscape/paths.local.json",
    "orfscape/tools/prodigal.exe",      # GPL-3 binary: install it yourself
    "orfscape/core/sjh.py",             # project-specific result pack, not part of the pipeline
    "scripts/_*.py",
    "base_view.png",
    "chain_after_e2.txt",
    "docs/GUI_CHANGELOG.md",
    "docs/GUI_BACKLOG.md",
]
# 发布包里不允许出现的内容（模式本身分段拼接，避免扫描器命中本文件）
FORBIDDEN = [
    (re.compile(r"[A-Za-z]:[\\/](?:JKY|Users|KL\.)", re.I), "personal absolute path"),
    (re.compile("JKY" + "_work|JKY" + "教务", re.I), "internal workspace path"),
    (re.compile("rsc" + "5385459"), "credential"),
    (re.compile(r"password\s*[:=]\s*['\"][^'\"]+['\"]", re.I), "hard-coded password"),
]
TEXT_EXT = {".py", ".md", ".txt", ".json", ".cff", ".cfg", ".ini", ".bat",
            ".yml", ".yaml", ".toml", ".tsv", ".csv"}


def _skip(path):
    for part in path.parts:
        if part in SKIP_DIRS:
            return True
    rel = path.relative_to(ROOT).as_posix()
    for pat in SKIP_GLOBS:
        if path.match(pat):
            return True
        if Path(rel).match(pat):
            return True
    return False


def collect():
    files = []
    for p in sorted(ROOT.rglob("*")):
        if p.is_dir() or _skip(p):
            continue
        files.append(p)
    return files


def scan(files):
    """返回 (问题列表)。只扫文本文件。"""
    problems = []
    me = Path(__file__).resolve()
    for p in files:
        if p.suffix.lower() not in TEXT_EXT or p.resolve() == me:
            continue
        try:
            text = p.read_text(encoding="utf-8", errors="replace")
        except OSError as exc:
            problems.append((p, "unreadable: %s" % exc))
            continue
        for rx, why in FORBIDDEN:
            m = rx.search(text)
            if m:
                line = text[:m.start()].count("\n") + 1
                problems.append((p, "%s at line %d: %s" % (why, line, m.group(0)[:60])))
    return problems


def version():
    init = ROOT / "orfscape" / "__init__.py"
    for line in init.read_text(encoding="utf-8").splitlines():
        if line.startswith("__version__"):
            return line.split("=", 1)[1].strip().strip('"\'')
    return "0.0.0"


def build(ver, files):
    DIST.mkdir(exist_ok=True)
    out = DIST / ("ORFscape_v%s.zip" % ver)
    arc_root = "ORFscape"
    with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED, compresslevel=9) as zf:
        for p in files:
            zf.write(p, arcname="%s/%s" % (arc_root, p.relative_to(ROOT).as_posix()))
    size = out.stat().st_size
    digest = hashlib.sha256(out.read_bytes()).hexdigest()
    checks = DIST / "SHA256SUMS.txt"
    checks.write_text("%s  %s\n" % (digest, out.name), encoding="utf-8")
    print("archive : %s" % out)
    print("files   : %d" % len(files))
    print("size    : %.1f MB" % (size / 1048576.0))
    print("sha256  : %s" % digest)
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description="ORFscape release builder")
    ap.add_argument("--check", action="store_true", help="scan only")
    args = ap.parse_args(argv)

    files = collect()
    problems = scan(files)
    if problems:
        print("REFUSING TO BUILD — %d problem(s):" % len(problems))
        for p, why in problems[:40]:
            print("  %s: %s" % (p.relative_to(ROOT), why))
        return 1
    ver = version()
    print("version : %s" % ver)
    print("files   : %d (clean)" % len(files))
    if args.check:
        return 0
    build(ver, files)
    return 0


if __name__ == "__main__":
    sys.exit(main())
