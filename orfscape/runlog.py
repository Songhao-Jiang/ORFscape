# -*- coding: utf-8 -*-
"""事件流 → 运行状态（GUI 的模型层，可在无显示环境下单测）。

    st = RunState()
    for rec in read_events(path):
        st.feed(rec)
    st.overall()      # 0..1 加权总进度
    st.eta_sec()      # 预计剩余秒数（无估计时返回 None）

设计要点
--------
- 进度按**阶段权重**加权，而不是按阶段个数，否则 E2（占大头）与 P8a 权重相同。
- ETA 优先用当前阶段实测速率（done/elapsed），不足两个采样点时退回
  ``run.start`` 给出的先验估计；再不足则用整体速率外推。
- 只读、纯函数式推进，任何字段缺失都不抛异常（日志可能被截断）。
"""
import datetime
import json
import os
import re

STATUS_PENDING = "pending"
STATUS_RUNNING = "running"
STATUS_OK = "ok"
STATUS_FAIL = "fail"
STATUS_SKIP = "skip"
STATUS_CHECK_FAIL = "check-fail"
STATUS_ABORTED = "aborted"

TERMINAL = {STATUS_OK, STATUS_FAIL, STATUS_SKIP, STATUS_CHECK_FAIL,
            STATUS_ABORTED}

# 阶段顺序（GUI 看板固定顺序；orchestrator 的 run.start 会覆盖）
DEFAULT_ORDER = ["P0", "P1", "D0", "C", "E1", "QC", "E2", "P6", "P6b", "P7",
                 "P8a", "P8b", "P8c"]


class Stage(object):
    """单个阶段的运行态。"""

    __slots__ = ("sid", "title_key", "weight", "est_sec", "status", "t0", "t1",
                 "last_t", "done", "total", "unit", "psm", "msg", "check",
                 "cmd", "err")

    def __init__(self, sid, title_key=None, weight=1.0, est_sec=None):
        self.sid = sid
        self.title_key = title_key or ("stage." + sid)
        self.weight = float(weight or 1.0)
        self.est_sec = est_sec
        self.status = STATUS_PENDING
        self.t0 = self.t1 = None
        self.last_t = None
        self.done = 0
        self.total = None
        self.unit = "task"
        self.psm = 0
        self.msg = None          # (key, args)
        self.check = None        # {"status":..., "fails":n, "warns":n}
        self.cmd = None
        self.err = None

    # ---------------------------------------------------------------- 派生值
    @property
    def elapsed(self):
        """已跑秒数。运行中用最后一条事件的时间（读旧日志也成立）。"""
        if self.t0 is None:
            return 0.0
        if self.t1 is not None:
            return self.t1 - self.t0
        return max(0.0, (self.last_t or _now()) - self.t0)

    def fraction(self):
        """阶段完成度 0..1。"""
        if self.status in (STATUS_OK, STATUS_SKIP):
            return 1.0
        if self.status in (STATUS_ABORTED, STATUS_FAIL, STATUS_CHECK_FAIL):
            # 失败/未跑的阶段没做完，进度条不能算满；有子进度就按比例
            return (max(0.0, min(1.0, float(self.done or 0) / float(self.total)))
                    if self.total else 0.0)
        if self.total:
            return max(0.0, min(1.0, float(self.done or 0) / float(self.total)))
        return 0.0

    def rate(self):
        """单位/秒（有子进度且已跑够 3 s 才有意义）。"""
        el = self.elapsed
        if not self.total or self.done <= 0 or el < 3:
            return None
        return self.done / el

    def remaining_sec(self):
        """阶段剩余秒数（可 None）。"""
        if self.status in TERMINAL:
            return 0.0
        if self.status != STATUS_RUNNING:
            return self.est_sec
        r = self.rate()
        if r and self.total:
            return max(0.0, (self.total - self.done) / r)
        if self.est_sec:
            f = self.fraction()
            return max(0.0, self.est_sec * (1.0 - f))
        return None

    def to_dict(self):
        return {"sid": self.sid, "status": self.status, "done": self.done,
                "total": self.total, "elapsed": round(self.elapsed, 1),
                "check": self.check, "psm": self.psm}


def _now():
    return datetime.datetime.now().timestamp()


class RunState(object):
    """整条流程的运行态。"""

    def __init__(self):
        self.order = list(DEFAULT_ORDER)
        self.stages = {}
        self.started = None
        self.ended = None
        self.status = STATUS_PENDING
        self.raw = None
        self.species = None
        self.desc = None
        self.workdir = None
        self.lang = None
        self.error = None
        self._eta_ema = None
        self.n_events = 0

    # ---------------------------------------------------------------- 构建
    def ensure(self, sid, **kw):
        st = self.stages.get(sid)
        if st is None:
            st = Stage(sid, **kw)
            self.stages[sid] = st
            if sid not in self.order:
                self.order.append(sid)
        return st

    def feed(self, rec):
        """吃一条事件；返回 True 表示状态有变化。"""
        if not isinstance(rec, dict):
            return False
        ev = rec.get("ev")
        if not ev:
            return False
        self.n_events += 1
        sid = rec.get("stage")
        t = rec.get("t") or _now()

        if ev == "run.start":
            self.started = t
            self.status = STATUS_RUNNING
            self.raw = rec.get("raw")
            self.species = rec.get("species")
            self.desc = rec.get("desc")
            self.workdir = rec.get("workdir")
            self.lang = rec.get("lang")
            stages = rec.get("stages") or []
            if stages:
                self.order = [s.get("id") for s in stages if s.get("id")] or self.order
                for s in stages:
                    st = self.ensure(s.get("id"), title_key=s.get("title_key"),
                                     weight=s.get("weight") or 1.0,
                                     est_sec=s.get("est_sec"))
                    st.weight = float(s.get("weight") or 1.0)
                    if s.get("est_sec"):
                        st.est_sec = float(s["est_sec"])
            return True

        if ev == "run.end":
            self.ended = t
            self.status = rec.get("status") or STATUS_OK
            self.error = rec.get("error")
            ok = self.status == STATUS_OK
            for st in self.stages.values():
                if st.status == STATUS_RUNNING:
                    st.status = STATUS_OK if ok else STATUS_FAIL
                    st.t1 = st.t1 or t
                elif st.status == STATUS_PENDING:
                    # 失败/中断时没跑到的阶段标 aborted，不能算作已完成
                    st.status = STATUS_SKIP if ok else STATUS_ABORTED
                    st.t1 = t
            return True

        if ev == "error":
            self.error = rec.get("text") or rec.get("line") or "error"
            return True

        if not sid:
            return True  # 无阶段上下文的事件（如全局提示）也算变化

        st = self.ensure(sid)
        if t and (st.last_t is None or t > st.last_t):
            # 只前进不后退：日志行可能带更早的时间戳（乱序/补写）
            st.last_t = t

        if ev == "stage.start":
            st.status = STATUS_RUNNING
            st.t0 = t
            st.t1 = None
            st.cmd = rec.get("cmd")
            if rec.get("est_sec"):
                st.est_sec = float(rec["est_sec"])
            return True

        if ev == "stage.enter":
            if st.t0 is None:
                st.t0 = t
                st.status = STATUS_RUNNING
            if rec.get("total"):
                st.total = int(rec["total"])
            return True

        if ev == "progress":
            if st.status == STATUS_PENDING:
                st.status = STATUS_RUNNING
                st.t0 = st.t0 or t
            if rec.get("total"):
                st.total = int(rec["total"])
            if rec.get("done") is not None:
                st.done = int(rec["done"])
            if rec.get("unit"):
                st.unit = rec["unit"]
            return True

        if ev == "metric":
            key = rec.get("key")
            if key == "psm":
                st.psm = int(rec.get("value") or 0)
            elif key == "total" and rec.get("value"):
                st.total = int(rec["value"])
            return True

        if ev == "msg":
            st.msg = (rec.get("key"), rec.get("args") or {})
            return True

        if ev == "stage.end":
            st.status = rec.get("status") or STATUS_OK
            # sec 与 (t1-t0) 可能差几毫秒，用上报值兜底：t0 反推
            if rec.get("sec") is not None:
                try:
                    st.t1 = t
                    st.t0 = t - float(rec["sec"])
                except (TypeError, ValueError):
                    st.t1 = t
            else:
                st.t1 = t
            if st.total and st.status == STATUS_OK:
                st.done = st.total
            return True

        if ev == "stage.check":
            st.check = {"status": rec.get("status"), "fails": rec.get("fails") or 0,
                        "warns": rec.get("warns") or 0}
            if (rec.get("status") or "").lower() in ("fail", "error"):
                st.status = STATUS_CHECK_FAIL
            return True

        return True

    # ---------------------------------------------------------------- 汇总
    def current(self):
        """正在跑的阶段（多个时取顺序最前的），否则最后一个已结束的。"""
        for sid in self.order:
            st = self.stages.get(sid)
            if st and st.status == STATUS_RUNNING:
                return st
        done = [self.stages[sid] for sid in self.order
                if sid in self.stages and self.stages[sid].status in TERMINAL]
        return done[-1] if done else None

    def failed(self):
        """出问题的阶段（fail / check-fail / stopped），没有返回 None。"""
        for sid in self.order:
            st = self.stages.get(sid)
            if st is not None and st.status in (STATUS_FAIL, STATUS_CHECK_FAIL,
                                                "stopped"):
                return st
        return None

    def overall(self):
        """加权总进度 0..1。"""
        if not self.order:
            return 0.0
        tot = 0.0
        acc = 0.0
        for sid in self.order:
            st = self.stages.get(sid)
            w = st.weight if st else 1.0
            tot += w
            if st is None:
                continue
            if st.status == STATUS_SKIP:
                acc += w                      # 跳过的阶段按已完成计
            else:
                acc += w * st.fraction()
        return max(0.0, min(1.0, acc / tot)) if tot else 0.0

    def eta_sec(self, smooth=True):
        """预计剩余秒数；完全无信息时返回 None。"""
        if self.status in TERMINAL or self.ended:
            return 0.0
        rem = 0.0
        unknown = False
        for sid in self.order:
            st = self.stages.get(sid)
            if st is None:
                unknown = True
                continue
            r = st.remaining_sec()
            if r is None:
                unknown = True
            else:
                rem += r
        if unknown and rem == 0.0 and self.started:
            f = self.overall()
            if f > 0.02:
                rem = max(0.0, (1.0 - f) / f * (_now() - self.started))
            else:
                return None
        if not smooth:
            return rem
        if self._eta_ema is None:
            self._eta_ema = rem
        else:
            # 只允许缓慢变化，避免数字乱跳
            self._eta_ema = 0.75 * self._eta_ema + 0.25 * rem
        return self._eta_ema

    def psm_total(self):
        """本轮到目前的 PSM 总数（各阶段累加）。"""
        return sum(int(getattr(st, "psm", 0) or 0) for st in self.stages.values())

    def counts(self):
        c = {"ok": 0, "running": 0, "pending": 0, "fail": 0, "skip": 0,
             "check-fail": 0, "aborted": 0, "stopped": 0}
        for sid in self.order:
            st = self.stages.get(sid)
            if st:
                c[st.status] = c.get(st.status, 0) + 1
        return c

    def rows(self):
        """给阶段表用的行（顺序固定）。"""
        out = []
        for sid in self.order:
            st = self.stages.get(sid) or Stage(sid)
            out.append(st)
        return out


# ---------------------------------------------------------------- 读文件
def read_events(path, from_offset=0):
    """从 offset 起读 JSONL，返回 (records, new_offset)。

    半行（进程正在写）不消费：offset 停在最后一条完整行之后。
    """
    out = []
    if not path or not os.path.isfile(path):
        return out, from_offset
    with open(path, "rb") as fh:
        fh.seek(from_offset, os.SEEK_SET)
        data = fh.read()
    if not data:
        return out, from_offset
    cut = data.rfind(b"\n")
    if cut < 0:
        return out, from_offset
    chunk, consumed = data[:cut + 1], cut + 1
    for raw in chunk.splitlines():
        raw = raw.strip()
        if not raw:
            continue
        try:
            rec = json.loads(raw.decode("utf-8", "replace"))
        except ValueError:
            continue
        if isinstance(rec, dict):
            out.append(rec)
    return out, from_offset + consumed


# ---------------------------------------------------------------- 兼容旧文本日志
_RE_STEP = re.compile(r"^>>>\s+(\w+):")
_RE_RC = re.compile(r"^\s*rc=(-?\d+)\s+([\d.]+)s")
_RE_SKIP = re.compile(r"^跳过\s+(\w+)")


class LegacyTail(object):
    """解析 run_pipeline.py 写出的 run_*.log（无 JSONL 时的降级路径）。

    只恢复"哪些阶段开始/结束/失败"，没有子进度，ETA 会用整体速率外推。
    """

    def __init__(self, state=None):
        self.state = state or RunState()
        self.lines = 0

    def feed_lines(self, lines):
        cur = None
        for line in lines:
            self.lines += 1
            m = _RE_STEP.match(line)
            if m:
                cur = m.group(1)
                st = self.state.ensure(cur)
                st.status = STATUS_RUNNING
                st.t0 = _now()
                continue
            if cur and _RE_RC.match(line):
                m = _RE_RC.match(line)
                st = self.state.stages.get(cur)
                if st:
                    st.t1 = _now()
                    st.status = STATUS_OK if int(m.group(1)) == 0 else STATUS_FAIL
                cur = None
        self.state.started = self.state.started or _now()
        self.state.status = STATUS_RUNNING
        return self.state
