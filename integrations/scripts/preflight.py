#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""ORFscape preflight: environment, dependencies, inputs, disk.

    python integrations/scripts/preflight.py --workdir orfscape_projects/<proj> --raw <raw folder>

Any FAIL means do not run the stages: ORFscape's failure mode is "it finished but
produced nothing", which costs hours.

Tool locations come from ORFSCAPE_<KEY> environment variables or from
`orfscape/paths.local.json` (copy `paths.local.example.json`).
"""
import argparse
import importlib.util
import json
import os
import shutil
import sys
from pathlib import Path

MIN_PY = (3, 9)


def _local_paths():
    """Read `orfscape/paths.local.json` from the package if present (optional)."""
    for base in Path(__file__).resolve().parents:
        for cand in (base / "orfscape" / "paths.local.json", base / "paths.local.json"):
            if cand.is_file():
                try:
                    with cand.open(encoding="utf-8") as fh:
                        data = json.load(fh)
                    if isinstance(data, dict):
                        return data
                except (OSError, ValueError):
                    return {}
    return {}


_LP = _local_paths()
# Priority: environment variable > paths.local.json > empty (a FAIL tells the user what to set)
DEFAULT_PFIND = os.path.join(str(_LP.get("pfind") or ""), "bin") if _LP.get("pfind") else ""
DEFAULT_ISOPSV = str(_LP.get("isopsv") or "")


def _p(level, msg):
    print("  [%s] %s" % (level, msg))


def _find_orfscape(start):
    """Walk up looking for a directory that contains orfscape/__init__.py."""
    cur = Path(start).resolve()
    for _ in range(6):
        for cand in (cur / "orfscape" / "__init__.py", cur / "ORFscape" / "orfscape" / "__init__.py"):
            if cand.is_file():
                return cand.parent.parent
        if cur.parent == cur:
            break
        cur = cur.parent
    return None


def check_python():
    ok = sys.version_info[:2] >= MIN_PY
    _p("PASS" if ok else "FAIL", "Python %d.%d (requires >=%d.%d)"
       % (sys.version_info[0], sys.version_info[1], MIN_PY[0], MIN_PY[1]))
    return ok


def check_orfscape(workdir):
    if importlib.util.find_spec("orfscape") is not None:
        _p("PASS", "orfscape is importable")
        return True
    root = _find_orfscape(workdir) or _find_orfscape(Path(__file__).parent)
    if root:
        _p("PASS", "orfscape found at %s (run `cd <that dir> && python -m orfscape`)" % root)
        return True
    _p("FAIL", "orfscape package not found: pip install -e ORFscape, or add the repo root to PYTHONPATH")
    return False


def check_pfind():
    path = os.environ.get("ORFSCAPE_PFIND") or DEFAULT_PFIND
    if not path or not os.path.isdir(path):
        _p("FAIL", "pFind directory missing: %s (set ORFSCAPE_PFIND)" % (path or "<unset>"))
        return False
    exes = [p for p in Path(path).glob("pFind*") if p.is_file()]
    _p("PASS" if exes else "WARN",
       "pFind directory %s (%d executables found)" % (path, len(exes)))
    return bool(exes)


def check_isopsv():
    path = os.environ.get("ORFSCAPE_ISOPSV") or DEFAULT_ISOPSV
    if not path or not os.path.isdir(path):
        _p("FAIL", "isoPSV directory missing: %s (set ORFSCAPE_ISOPSV)" % (path or "<unset>"))
        return False
    _p("PASS", "isoPSV directory present: %s" % path)
    return True


def check_workdir(workdir):
    wd = Path(workdir).resolve()
    if not wd.is_dir():
        _p("FAIL", "project directory does not exist: %s" % wd)
        return False
    probe = wd / ".write_probe"
    try:
        probe.write_text("ok", encoding="utf-8")
        probe.unlink()
    except OSError as exc:
        _p("FAIL", "project directory is not writable: %s" % exc)
        return False
    free_gb = shutil.disk_usage(str(wd)).free / (1024 ** 3)
    _p("PASS" if free_gb >= 20 else "WARN",
       "project directory writable, %.1f GB free (>=20 GB recommended)" % free_gb)
    return True


def check_raw(raw):
    if not raw:
        _p("WARN", "no --raw given, skipping raw data check")
        return True
    root = Path(raw).resolve()
    if not root.is_dir():
        _p("FAIL", "raw directory does not exist: %s" % root)
        return False
    n = 0
    for ext in ("*.raw", "*.d", "*.wiff", "*.dta"):
        n += sum(1 for _ in root.rglob(ext))
    _p("PASS" if n else "FAIL", "raw directory %s: %d raw file(s)" % (root, n))
    return n > 0


def check_genome(workdir):
    base = Path(workdir).resolve() / "01_genome"
    if not base.is_dir():
        _p("WARN", "01_genome not present yet (normal before P1)")
        return True
    species = [d for d in base.iterdir() if d.is_dir()]
    if not species:
        _p("WARN", "no species directory under 01_genome")
        return True
    ok = True
    for sp in species:
        miss = [f for f in ("genomic.fna", "genomic.gff", "protein.faa")
                if not (sp / f).is_file()]
        if miss:
            _p("FAIL", "%s missing: %s" % (sp.name, ", ".join(miss)))
            ok = False
        else:
            _p("PASS", "%s genome triplet present" % sp.name)
    return ok


def main(argv=None):
    ap = argparse.ArgumentParser(description="ORFscape preflight")
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--raw", default=None)
    args = ap.parse_args(argv)

    print("ORFscape preflight")
    results = [
        check_python(),
        check_orfscape(args.workdir),
        check_pfind(),
        check_isopsv(),
        check_workdir(args.workdir),
        check_raw(args.raw),
        check_genome(args.workdir),
    ]
    bad = results.count(False)
    print("\nresult: %s (%d failure(s))" % ("FAIL" if bad else "OK", bad))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
