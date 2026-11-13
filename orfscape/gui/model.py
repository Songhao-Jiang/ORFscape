# -*- coding: utf-8 -*-
"""GUI 的输入模型与偏好持久化（无 tkinter 依赖，可单测）。

偏好只写一个固定位置：``~/.orfscape_gui.json``。路径在导入时算好、规范化，
不接受调用方传路径；读写都走 pathlib，且写前确认目标仍在用户主目录。
"""
import json
import os
from pathlib import Path

_REAL_HOME = Path(os.path.realpath(os.path.expanduser("~")))
SETTINGS = _REAL_HOME / ".orfscape_gui.json"
# 改名（2026-11-13）前叫 .prokpg_gui.json：新文件不存在时把旧偏好迁过来一次
_LEGACY = _REAL_HOME / (".prokpg" + "_gui.json")

DEFAULTS = {
    "raws": [],
    "species": "",
    "desc": "",
    "workdir": "",
    "threads": 16,
    "workers": 4,
    "dbs": "aug,6f",
    "skip_p6b": False,
    "skip_qc": False,
    "force": False,
    "lang": "en",
    "advanced": False,
    "geometry": "",
    "theme": "light",
    "confirm": True,
    "font_scale": 1.0,
    "filter": "",
    "log_level": 0,
}


def _settings_file():
    """规范化后的偏好文件；必须正好位于用户主目录下。"""
    p = Path(os.path.realpath(str(SETTINGS)))
    if p.parent != _REAL_HOME or p.suffix != ".json":
        return None
    return p


def _migrate_legacy():
    """旧偏好文件存在而新文件不存在时，复制一份（只做一次）。"""
    p = _settings_file()
    if p is None or p.is_file():
        return
    try:
        if _LEGACY.is_file() and _LEGACY.parent == _REAL_HOME:
            p.write_bytes(_LEGACY.read_bytes())
    except OSError:
        pass


def load_settings():
    data = dict(DEFAULTS)
    _migrate_legacy()
    p = _settings_file()
    if p is None or not p.is_file():
        return data
    try:
        saved = json.loads(p.read_text(encoding="utf-8"))
        if isinstance(saved, dict):
            for k in DEFAULTS:
                if k in saved and isinstance(saved[k], type(DEFAULTS[k])):
                    data[k] = saved[k]
            if isinstance(saved.get("raws"), list):
                data["raws"] = [str(x) for x in saved["raws"] if isinstance(x, str)]
    except (OSError, ValueError):
        pass
    return data


def save_settings(data):
    """尽力而为地保存（写不进去不影响使用）。"""
    p = _settings_file()
    if p is None:
        return False
    out = {k: data.get(k, DEFAULTS[k]) for k in DEFAULTS}
    try:
        p.write_text(json.dumps(out, ensure_ascii=False, indent=1), encoding="utf-8")
        return True
    except OSError:
        return False


def common_raw_root(files):
    """多个 raw 的最近公共目录；单文件返回其父目录。"""
    files = [os.path.abspath(f) for f in files if f]
    if not files:
        return ""
    dirs = [f if os.path.isdir(f) else os.path.dirname(f) for f in files]
    if len(dirs) == 1:
        return dirs[0]
    try:
        return os.path.commonpath(dirs)
    except ValueError:
        return ""


def raw_stats(files):
    n = 0
    total = 0
    for f in files or []:
        if os.path.isdir(f):
            for root, _dirs, names in os.walk(f):
                for name in names:
                    if name.lower().endswith(".raw"):
                        n += 1
                        try:
                            total += os.path.getsize(os.path.join(root, name))
                        except OSError:
                            pass
        else:
            # 不存在也算一条（用户刚选完、盘符还没挂上时不该显示"未选择"）
            n += 1
            try:
                total += os.path.getsize(f)
            except OSError:
                pass
    return n, total


def build_run_args(form, workdir, lang):
    """把表单转成 ``orfscape run`` 的参数表（不含解释器与 ``-m orfscape``）。"""
    args = ["run", "--workdir", workdir]
    root = common_raw_root(form.get("raws") or [])
    if root:
        args += ["--raw", root]
    if form.get("species"):
        args += ["--species-name", str(form["species"]).strip()]
    if form.get("desc"):
        args += ["--desc", str(form["desc"])]
    if form.get("threads"):
        args += ["--thread", str(int(form["threads"]))]
    if form.get("workers"):
        args += ["--workers", str(int(form["workers"]))]
    if form.get("dbs"):
        args += ["--dbs", str(form["dbs"])]
    skip = []
    if form.get("skip_p6b"):
        skip.append("P6b")
    if form.get("skip_qc"):
        skip.append("QC")
    if skip:
        args += ["--skip", ",".join(skip)]
    if form.get("force"):
        args += ["--force"]
    args += ["--lang", lang]
    return args


def line_level(line):
    """日志行级别：0=普通 1=警告 2=错误。"""
    low = (line or "").lower()
    if ("fail" in low or "error" in low or "traceback" in low
            or "✕" in line or "失败" in line):
        return 2
    if "warn" in low or "⚠" in line or "警告" in line:
        return 1
    return 0


def filter_log(entries, pattern="", min_level=0, stage_ids=()):
    """按阶段/关键词/级别过滤日志。

    entries: [(line, stage)]
    pattern: 正好是阶段号时按阶段**精确**匹配（P6 不会带出 P6b），
             否则按子串匹配。
    返回 [(line, level)]。
    """
    pat = (pattern or "").strip().lower()
    by_stage = pat in {s.lower() for s in stage_ids}
    out = []
    for line, stg in entries:
        lv = line_level(line)
        if min_level and lv < min_level:
            continue
        if by_stage:
            if (stg or "").lower() != pat:
                continue
        elif pat and pat not in (line or "").lower():
            continue
        out.append((line, lv))
    return out


def resume_stage(state):
    """续跑起点：第一个不是 ok/skip 的阶段；全完成返回 None。"""
    if state is None:
        return None
    for sid in state.order:
        st = state.stages.get(sid)
        if st is None:
            return sid
        if st.status not in ("ok", "skip"):
            return sid
    return None


def stage_tooltip(stage, est_sec=None):
    """阶段悬浮提示的内容（结构化，交给界面翻译）。返回 [(key, args)]。"""
    from ..estimator import fmt_sec
    out = []
    if stage is None:
        return out
    out.append(("tip.status", {"v": stage.status}))
    if stage.total:
        out.append(("tip.progress", {"v": "%d/%d %s" % (stage.done or 0,
                                                        stage.total, stage.unit)}))
    elif stage.psm:
        out.append(("tip.progress", {"v": "PSM {:,}".format(stage.psm)}))
    if stage.elapsed:
        out.append(("tip.elapsed", {"v": fmt_sec(stage.elapsed),
                                    "e": fmt_sec(est_sec) if est_sec else "—"}))
    if stage.check:
        out.append(("tip.check", {"f": stage.check.get("fails") or 0,
                                  "w": stage.check.get("warns") or 0}))
    return out


def latest_log(workdir):
    """该项目最新的一份 `run_*.log`（没有返回 None）。"""
    if not workdir:
        return None
    logs = Path(os.path.realpath(workdir)) / "00_logs"
    if not logs.is_dir():
        return None
    try:
        cands = sorted((p for p in logs.glob("run_*.log") if p.is_file()),
                       key=lambda p: p.stat().st_mtime)
    except OSError:
        return None
    return cands[-1] if cands else None


def find_occurrences(text, needle):
    """大小写不敏感地找全部命中，返回 [(start, end)]（用于计数与单测）。"""
    text = text or ""
    needle = needle or ""
    if not needle:
        return []
    low, pat = text.lower(), needle.lower()
    out, pos = [], low.find(pat)
    while pos >= 0 and len(out) < 5000:
        out.append((pos, pos + len(pat)))
        pos = low.find(pat, pos + len(pat))
    return out


def summary_rows(state, est=None):
    """把运行态摊平成表格行（导出与单测共用）。

    列：stage / status / elapsed_s / est_s / check_fails / check_warns / psm
    """
    rows = []
    for sid in (state.order if state else []):
        st = state.stages.get(sid)
        if st is None:
            rows.append([sid, "pending", "", "", "", "", ""])
            continue
        est_s = ((est or {}).get("est_sec", {}) or {}).get(sid) or {}
        check = st.check or {}
        rows.append([
            sid,
            st.status,
            "%.1f" % st.elapsed if st.elapsed else "",
            ("%.0f" % est_s["sec"]) if est_s.get("sec") else "",
            check.get("fails", ""),
            check.get("warns", ""),
            st.psm or "",
        ])
    return rows


def export_summary(state, workdir, stamp=None, est=None):
    """写 `00_logs/summary_<stamp>.csv`，返回 Path；失败抛 OSError。"""
    import csv
    import io
    import time
    stamp = stamp or time.strftime("%Y%m%d_%H%M%S")
    base = Path(os.path.realpath(workdir))
    logs = base / "00_logs"
    logs.mkdir(parents=True, exist_ok=True)
    target = logs / ("summary_%s.csv" % stamp)
    if target.parent != logs:
        raise ValueError("unsafe summary path")
    buf = io.StringIO(newline="")
    w = csv.writer(buf)
    w.writerow(["project", state.workdir or str(base)])
    w.writerow(["started", time.strftime("%Y-%m-%d %H:%M:%S",
                                         time.localtime(state.started))
                if state.started else ""])
    w.writerow(["finished", time.strftime("%Y-%m-%d %H:%M:%S",
                                          time.localtime(state.ended))
                if state.ended else ""])
    w.writerow(["run_status", state.status])
    w.writerow(["overall", "%.3f" % state.overall()])
    w.writerow(["psm_total", state.psm_total()])
    w.writerow([])
    w.writerow(["stage", "status", "elapsed_s", "est_s", "check_fails",
                "check_warns", "psm"])
    for row in summary_rows(state, est):
        w.writerow(row)
    target.write_text(buf.getvalue(), encoding="utf-8-sig")
    return target
