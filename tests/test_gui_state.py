# -*- coding: utf-8 -*-
"""事件状态机的随机推进测试（无窗口）。

    python ORFscape/tests/test_gui_state.py [--rounds 200]

对随机生成的事件序列断言不变量：
- 0 ≤ 总进度 ≤ 1，且不会越过 1；
- 只要没有阶段被"重跑"（状态只前进），总进度不下降；
- 只有一个阶段处于 running 时 `current()` 必须返回它；
- 正常结束 → 1.0 且 ETA=0；失败结束 → <1.0 且未跑阶段是 aborted；
- ETA 恒 ≥0 或 None；PSM 汇总等于各阶段 metric 之和；
- 任意顺序喂事件都不抛异常。
"""
import argparse
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

from orfscape import runlog                                  # noqa: E402

ORDER = ["P0", "P1", "D0", "C", "E1", "QC", "E2", "P6", "P6b", "P7",
         "P8a", "P8b", "P8c"]


class Rng(object):
    """确定性小 PRNG：测试只要可复现，不需要安全随机。"""

    _M = (1 << 64) - 1
    _A = 6364136223846793005
    _C = 1442695040888963407

    def __init__(self, seed):
        self.s = (int(seed) * self._A + self._C) & self._M

    def next(self):
        self.s = (self.s * self._A + self._C) & self._M
        return (self.s >> 11) / float(1 << 53)

    def random(self):
        return self.next()

    def uniform(self, a, b):
        return a + (b - a) * self.next()

    def choice(self, seq):
        return seq[int(self.next() * len(seq)) % len(seq)]

    def sample(self, seq, k):
        pool, out = list(seq), []
        for _ in range(min(k, len(pool))):
            out.append(pool.pop(int(self.next() * len(pool)) % len(pool)))
        return out


def start_event(order, t0=0.0):
    return {"ev": "run.start", "t": t0, "workdir": "/x",
            "stages": [{"id": s, "weight": 1.0 + i, "est_sec": 60.0 + 10 * i}
                       for i, s in enumerate(order)]}


def gen_events(rng, order, outcome):
    """产生一段合法事件流：顺序推进，随机跳过/失败/无子进度。"""
    recs = [start_event(order)]
    t = 100.0
    fail_at = rng.choice(order) if outcome == "fail" else None
    for sid in order:
        if fail_at and sid == fail_at:
            recs.append({"ev": "stage.start", "stage": sid, "t": t})
            t += rng.uniform(1, 30)
            recs.append({"ev": "stage.end", "stage": sid, "status": "fail",
                         "sec": t - 100.0})
            recs.append({"ev": "run.end", "status": "fail", "t": t})
            return recs
        roll = rng.random()
        if roll < 0.15:
            recs.append({"ev": "stage.start", "stage": sid, "skipped": True})
            recs.append({"ev": "stage.end", "stage": sid, "status": "skip",
                         "sec": 0, "t": t})
            continue
        recs.append({"ev": "stage.start", "stage": sid, "t": t})
        total = rng.choice([0, 5, 40, 144])
        if total:
            steps = rng.sample(range(total + 1), k=min(4, total + 1))
            for d in sorted(steps):
                t += rng.uniform(0.5, 8)
                recs.append({"ev": "progress", "stage": sid, "done": d,
                             "total": total, "t": t})
        t += rng.uniform(1, 60)
        recs.append({"ev": "stage.end", "stage": sid, "status": "ok",
                     "sec": 10.0, "t": t})
        recs.append({"ev": "stage.check", "stage": sid, "status": "ok",
                     "fails": 0, "warns": rng.choice([0, 0, 1]), "t": t})
    recs.append({"ev": "run.end", "status": "ok", "t": t})
    return recs


def run_round(seed, outcome):
    rng = Rng(seed)
    order = ORDER[:]
    recs = gen_events(rng, order, outcome)
    st = runlog.RunState()
    prev = 0.0
    for i, rec in enumerate(recs):
        st.feed(rec)
        frac = st.overall()
        assert 0.0 <= frac <= 1.0 + 1e-9, "overall out of range: %r" % frac
        if outcome == "ok":
            assert frac >= prev - 1e-9, "progress went backwards at %d" % i
        prev = frac
        eta = st.eta_sec(smooth=False)
        assert eta is None or eta >= 0, "negative ETA: %r" % eta
        cur = st.current()
        running = [s for s in st.rows() if s.status == "running"]
        if running:
            assert cur is not None and cur.sid == running[0].sid, \
                "current() inconsistent"
        assert cur is None or cur.sid in [s.sid for s in st.rows()]
    if outcome == "ok":
        assert abs(st.overall() - 1.0) < 1e-9, "ok run should reach 1.0"
        assert st.eta_sec() == 0.0
        assert not [s for s in st.rows() if s.status in ("aborted", "pending")], \
            "ok run left unfinished stages"
    else:
        assert st.overall() < 1.0, "failed run should not reach 1.0"
        assert st.failed() is not None, "failed() should find the stage"
        assert any(s.status == "aborted" for s in st.rows()) or \
            st.failed().sid == st.order[-1], "unvisited stages should be aborted"


def test_garbage():
    """乱序/缺字段事件也不能抛异常。"""
    rng = Rng(7)
    st = runlog.RunState()
    for _ in range(500):
        rec = {"ev": rng.choice(["progress", "metric", "log", "stage.end",
                                 "stage.check", "msg", "run.end", "nope", None]),
               "stage": rng.choice(ORDER + [None, "ZZ"]),
               "done": rng.choice([None, 0, 3, -1]),
               "total": rng.choice([None, 0, 10]),
               "status": rng.choice([None, "ok", "fail", "skip"])}
        st.feed(rec)
        assert 0.0 <= st.overall() <= 1.0
    st.feed("not a dict")
    assert st.overall() >= 0.0


def test_psm_sum():
    st = runlog.RunState()
    st.feed(start_event(ORDER[:3]))
    st.feed({"ev": "metric", "stage": "E2", "key": "psm", "value": 120})
    st.feed({"ev": "metric", "stage": "P6", "key": "psm", "value": 30})
    assert st.psm_total() == 150
    st.feed({"ev": "metric", "stage": "P6", "key": "psm", "value": 5})
    assert st.psm_total() == 125, "psm metric is absolute per stage"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rounds", type=int, default=120)
    args = ap.parse_args()
    fails = 0
    for i in range(args.rounds):
        for outcome in ("ok", "fail"):
            try:
                run_round(i * 2 + (outcome == "fail"), outcome)
            except AssertionError as exc:
                fails += 1
                print("FAIL seed=%d outcome=%s: %s" % (i, outcome, exc))
    test_garbage()
    test_psm_sum()
    print("%d rounds x2 outcomes, %d failure(s)" % (args.rounds, fails))
    return 0 if not fails else 1


if __name__ == "__main__":
    sys.exit(main())
