# -*- coding: utf-8 -*-
"""结构化进度事件：CLI/agent 写，GUI 读。

设计
----
- 事件以 JSONL（每行一个对象）追加到 ``$ORFSCAPE_EVENTS`` 指向的文件；
  变量未设置时全部 no-op，对命令行直跑零影响。
- 每个事件固定带 ``t``（epoch 秒）、``ts``（本地 ISO）、``pid``、``ev``；
  有阶段上下文时再带 ``stage``。
- 文案不写死：需要人读的提示用 ``msg(key, **args)``，GUI 用自己的
  i18n 目录翻译，日志里只存 key 与参数。
- 同一进程内所有线程共用一个句柄（P6b 是多线程），写操作加锁。

事件类型
--------
run.start / stage.enter / stage.start / progress / metric / msg / log /
stage.end / stage.check / run.end / error

GUI 侧解析见 :mod:`orfscape.runlog`。
"""
import json
import os
import re
import threading
import time

ENV_PATH = "ORFSCAPE_EVENTS"

_LOCK = threading.Lock()
_HANDLE = None
_HANDLE_PATH = None
_STAGE = None
_TOTAL = None
_LAST_TICK = (None, 0.0)   # (done, t) —— 去重/限流
_PSM = 0

# 进度嗅探：控制台既有输出里的 "[12/36]" 模式
_RE_TICK = re.compile(r"\[\s*(\d+)\s*/\s*(\d+)\s*\]")
_RE_TOTAL = re.compile(r"(?:任务数|tasks?)\s*[:：]\s*(\d+)")
_RE_PSM = re.compile(r"PSM\s*=\s*(\d+)")


def path():
    """当前事件文件路径（未启用返回 None）。"""
    p = os.environ.get(ENV_PATH, "").strip()
    return p or None


def enabled():
    return path() is not None


def _now():
    t = time.time()
    return t, time.strftime("%Y-%m-%dT%H:%M:%S", time.localtime(t))


def _fh():
    """懒打开事件文件；路径变化（父进程换项目）时重开。"""
    global _HANDLE, _HANDLE_PATH
    p = path()
    if not p:
        return None
    if _HANDLE is not None and _HANDLE_PATH == p:
        return _HANDLE
    if _HANDLE is not None:
        try:
            _HANDLE.close()
        except OSError:
            pass
    parent = os.path.dirname(os.path.abspath(p))
    if parent and not os.path.isdir(parent):
        os.makedirs(parent, exist_ok=True)
    _HANDLE = open(p, "a", encoding="utf-8", newline="\n")
    _HANDLE_PATH = p
    return _HANDLE


def emit(ev, **kw):
    """写一个事件；仅在 ORFSCAPE_EVENTS 生效时真正落盘。"""
    fh = _fh()
    if fh is None:
        return
    t, ts = _now()
    rec = {"t": t, "ts": ts, "pid": os.getpid(), "ev": ev}
    if _STAGE:
        rec["stage"] = _STAGE
    rec.update(kw)
    line = json.dumps(rec, ensure_ascii=False, separators=(",", ":"), default=str)
    with _LOCK:
        try:
            fh.write(line + "\n")
            fh.flush()
        except (OSError, ValueError):
            pass


def set_stage(stage, total=None):
    """设置当前阶段（子命令入口调用），并可预置任务总数。"""
    global _STAGE, _TOTAL, _LAST_TICK, _PSM
    _STAGE, _TOTAL, _LAST_TICK, _PSM = stage, total, (None, 0.0), 0
    if stage:
        emit("stage.enter", id=stage, total=total)


def tick(done, total=None, unit="task", **kw):
    """报一次子进度。done 相同且 1 s 内的重复调用直接丢弃。"""
    global _TOTAL, _LAST_TICK
    if total:
        _TOTAL = int(total)
    total = _TOTAL
    done = int(done)
    t = time.time()
    if _LAST_TICK[0] == done and (t - _LAST_TICK[1]) < 1.0:
        return
    _LAST_TICK = (done, t)
    emit("progress", done=done, total=total, unit=unit, **kw)


def metric(key, value, total=None, **kw):
    emit("metric", key=key, value=value, total=total, **kw)


def msg(key, **args):
    """给人看的提示：只存 key 与参数，由 GUI 翻译。"""
    emit("msg", key=key, args=args)


def log(line, **kw):
    """控制台一行的镜像（同时嗅探 [i/N] / PSM= / 任务数）。"""
    global _TOTAL, _PSM
    m = _RE_TICK.search(line)
    if m:
        tick(int(m.group(1)), int(m.group(2)))
    else:
        m = _RE_TOTAL.search(line)
        if m and _TOTAL is None:
            _TOTAL = int(m.group(1))
            emit("metric", key="total", value=_TOTAL, total=_TOTAL)
    m = _RE_PSM.search(line)
    if m:
        _PSM += int(m.group(1))
        metric("psm", _PSM)
    emit("log", line=line.rstrip("\n")[:4000], **kw)


def bytes_done(done, total, unit="B", **kw):
    """下载/拷贝类进度（字节）。"""
    emit("progress", done=int(done), total=int(total) if total else None,
         unit=unit, **kw)


def stage_end(status="ok", **kw):
    global _STAGE
    emit("stage.end", id=_STAGE, status=status, **kw)
    _STAGE = None


def close():
    global _HANDLE
    with _LOCK:
        if _HANDLE is not None:
            try:
                _HANDLE.close()
            except OSError:
                pass
        _HANDLE = None
