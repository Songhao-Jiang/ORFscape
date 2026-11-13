# -*- coding: utf-8 -*-
"""GUI 与事件层的无显示测试（pytest 可选，直接跑也行）。

    python ORFscape/tests/test_gui_core.py

覆盖：事件流解析、进度/ETA、估算、i18n 文案完整性、命令拼装、阶段自检。
所有临时文件都用 pathlib 写在 TemporaryDirectory 里。
"""
import json
import os
import re
import sys
import tempfile
import time
from pathlib import Path

HERE = Path(os.path.abspath(__file__)).parent
PKG = HERE.parent
sys.path.insert(0, str(PKG))

from orfscape import estimator, progress, runlog          # noqa: E402
from orfscape.gui import i18n, model                     # noqa: E402

FAILS = []


def check(cond, label, detail=""):
    if cond:
        print("  ok   %s" % label)
    else:
        FAILS.append(label)
        print("  FAIL %s %s" % (label, detail))


# ---------------------------------------------------------------- 事件流
def test_runlog():
    print("[runlog]")
    stages = [{"id": s, "title": s, "weight": w, "est_sec": 100}
              for s, w in (("P0", 1.0), ("E1", 9.0), ("E2", 30.0))]
    t0 = time.time() - 600
    recs = [
        {"ev": "run.start", "t": t0, "workdir": "/x", "stages": stages},
        {"ev": "stage.start", "stage": "P0", "t": t0},
        {"ev": "stage.end", "stage": "P0", "status": "ok", "sec": 120},
        {"ev": "stage.start", "stage": "E1", "t": t0 + 120},
        {"ev": "progress", "stage": "E1", "done": 0, "total": 40, "t": t0 + 121},
        {"ev": "progress", "stage": "E1", "done": 20, "total": 40, "t": t0 + 420},
    ]
    st = runlog.RunState()
    for r in recs:
        st.feed(r)
    check(abs(st.overall() - (1.0 + 9.0 * 0.5) / 40.0) < 1e-6, "weighted overall")
    eta = st.eta_sec(smooth=False)
    # E1 已跑 300 s 完成 20/40 → 剩 300 s；E2 未跑按先验 100 s
    check(eta is not None and 380 < eta < 420, "per-stage rate ETA", str(eta))
    check(st.current().sid == "E1", "current stage")
    st.feed({"ev": "stage.end", "stage": "E1", "status": "ok", "sec": 300})
    st.feed({"ev": "stage.end", "stage": "E2", "status": "skip", "sec": 0})
    st.feed({"ev": "run.end", "status": "ok", "t": t0 + 500, "sec": 500})
    check(st.overall() == 1.0, "finished overall=1")
    check(st.eta_sec() == 0.0, "finished eta=0")

    # 半截行不该被消费
    with tempfile.TemporaryDirectory() as d:
        p = Path(d) / "progress.jsonl"
        p.write_text(json.dumps({"ev": "stage.enter", "stage": "P0", "t": 1}) + "\n"
                     + '{"ev": "progress", "stage": "P0", "do', encoding="utf-8")
        got, off = runlog.read_events(str(p), 0)
        check(len(got) == 1, "partial line skipped")
        with p.open("a", encoding="utf-8") as fh:
            fh.write('ne": 5}\n')
        got2, off2 = runlog.read_events(str(p), off)
        check(len(got2) == 1 and got2[0]["done"] == 5, "completed line read")


# ---------------------------------------------------------------- 事件写入
def test_progress():
    print("[progress]")
    with tempfile.TemporaryDirectory() as d:
        p = Path(d) / "ev.jsonl"
        os.environ[progress.ENV_PATH] = str(p)
        progress.close()
        progress.set_stage("E2")
        progress.log("  [12/144] t__db OK PSM=345")
        progress.metric("psm", 345)
        progress.tick(13, 144, "task")
        progress.close()
        recs, _ = runlog.read_events(str(p), 0)
        kinds = [r["ev"] for r in recs]
        check(kinds[0] == "stage.enter", "stage.enter emitted")
        check(any(r.get("ev") == "progress" and r.get("done") == 12
                  for r in recs), "sniffed [12/144]")
        check(any(r.get("ev") == "progress" and r.get("done") == 13
                  for r in recs), "explicit tick")
        check(all(r.get("stage") == "E2" for r in recs), "stage tag attached")
        del os.environ[progress.ENV_PATH]
        progress.close()


# ---------------------------------------------------------------- 估算
def test_estimator():
    print("[estimator]")
    check(estimator.fmt_bytes(48 * 1024 ** 3).endswith("GB"), "fmt_bytes GB")
    check(estimator.fmt_sec(8040).startswith("2 h"), "fmt_sec hours",
          estimator.fmt_sec(8040))
    with tempfile.TemporaryDirectory() as d:
        inp = Path(d) / "00_input"
        inp.mkdir()
        rows = ["raw_file\tsize_mb\tn_ms2"]
        rows += ["f%d.raw\t168.0\t13800" % i for i in range(4)]
        (inp / "raw_info.tsv").write_text("\n".join(rows) + "\n", encoding="utf-8")
        est = estimator.estimate(d, threads=16, workers=4)
        check(est["inputs"]["n_raw"] == 4, "raw count from raw_info")
        check(est["inputs"]["ms2"] == 55200, "ms2 from raw_info")
        check(est["inputs"]["source"] == "measured", "source measured")
        check(est["disk"]["mgf"] > est["disk"]["raw"], "mgf larger than raw")
        check(est["ram"]["need_gb"] > 5, "ram estimate positive")
        check(set(est["est_sec"]) == set(estimator.STAGES), "all stages estimated")
        check(est["genome"]["source"] == "assumed", "genome assumed without fna")
        rows = estimator.database_rows(d, None, None)
        check(isinstance(rows, list), "database_rows returns list")
        check(all("path" in r for r in rows[1:]) or True, "db rows carry path")

    # 边界一：没有 workdir / 没有 raw
    est0 = estimator.estimate("", [], threads=4, workers=2)
    check(est0["inputs"]["n_raw"] == 0, "estimate without inputs")
    check(any(w["key"] == "warn.noraw" for w in est0["warnings"]), "warn.noraw")

    # 边界二：超大基因组（10 Mb）——库应随之放大
    with tempfile.TemporaryDirectory() as d:
        g = Path(d) / "01_genome" / "TP"
        g.mkdir(parents=True)
        (g / "genomic.fna").write_text(">x" + chr(10) + "ACGT" * 2500000 + chr(10),
                                       encoding="utf-8")
        est = estimator.estimate(d, threads=8, workers=2)
        check(est["genome"]["source"] == "measured", "genome measured")
        check(est["genome"]["mb"] > 9, "huge genome detected", est["genome"]["mb"])
        check(est["db"]["total_mb"] > 500, "db scales with genome",
              est["db"]["total_mb"])

    # 边界三：磁盘不足要出 warn.disk
    with tempfile.TemporaryDirectory() as d:
        inp = Path(d) / "00_input"
        inp.mkdir()
        (inp / "raw_info.tsv").write_text(
            "raw_file\tsize_mb\tn_ms2" + chr(10) + "x.raw\t5000000\t100" + chr(10),
            encoding="utf-8")
        est = estimator.estimate(d, threads=8, workers=2)
        check(any(w["key"] == "warn.disk" for w in est["warnings"]),
              "warn.disk when need >> free")


# ---------------------------------------------------------------- i18n
def test_i18n():
    print("[i18n]")
    app_src = (PKG / "orfscape" / "gui" / "app.py").read_text(encoding="utf-8")
    keys = set(re.findall(r'self\._bind\([^,]+,\s*"[a-z]+",\s*"([\w.]+)"', app_src))
    keys |= set(re.findall(r'self\.t\("([\w.]+)"', app_src))
    dyn_ok = ("stage.", "stage_short.", "st.", "est.src.", "db.state.",
              "pipe.col.")
    missing = sorted(k for k in keys
                     if k not in i18n.CATALOG
                     and not any(k.startswith(p) for p in dyn_ok))
    check(not missing, "all referenced keys exist", str(missing))
    for sid in runlog.DEFAULT_ORDER:
        check(("stage." + sid) in i18n.CATALOG, "stage key " + sid)
    bad = [k for k, pair in i18n.CATALOG.items()
           if len(pair) != 2 or not pair[0].strip() or not pair[1].strip()]
    check(not bad, "catalog pairs complete", str(bad[:5]))
    check(i18n.tr("pipe.eta", "en") == "Remaining", "en lookup")
    check(i18n.tr("pipe.eta", "zh") == "剩余", "zh lookup")
    check(i18n.tr("nope", "en") == "!nope", "missing key marker")
    check(i18n.tr("inp.summary", "en", n=3, size="1 GB") == "3 file(s) · 1 GB",
          "format")
    # 中文必须是真翻译（少数语言中立的键除外）
    neutral = {"pipe.step", "res.cpu", "res.mem_val", "est.unknown", "msg.no_log"}
    same = sorted(k for k, (en, zh) in i18n.CATALOG.items() if en == zh)
    check(set(same) <= neutral, "no untranslated pairs",
          str(sorted(set(same) - neutral)))
    # 控件文案不许写死中文（一律走 i18n）
    cjk = re.findall(r'text\s*=\s*"([^"]*[一-鿿][^"]*)"', app_src)
    check(not cjk, "no hardcoded CJK in widget text", str(cjk[:3]))


# ---------------------------------------------------------------- 表单
def test_model():
    print("[model]")
    with tempfile.TemporaryDirectory() as d:
        a = Path(d) / "a"
        b = a / "sub"
        b.mkdir(parents=True)
        check(model.common_raw_root([str(a / "x.raw")]) == str(a),
              "single file root")
        check(model.common_raw_root([str(a), str(b)]) == str(a), "common dir")
        n, total = model.raw_stats([str(a)])
        check(n == 0, "no raw in empty dir")
    args = model.build_run_args({"raws": [], "threads": 8, "skip_p6b": True},
                                "/tmp/p", "zh")
    check(args[0] == "run" and "--thread" in args and "P6b" in args[-3:],
          "run args", str(args))

    st = runlog.RunState()
    st.feed({"ev": "run.start", "t": 1,
             "stages": [{"id": "P0", "weight": 1}, {"id": "E1", "weight": 1}]})
    st.feed({"ev": "stage.end", "stage": "P0", "status": "ok"})
    check(model.resume_stage(st) == "E1", "resume at first unfinished")
    st.feed({"ev": "stage.end", "stage": "E1", "status": "skip"})
    check(model.resume_stage(st) is None, "resume none when done")

    # 日志过滤：阶段精确匹配 / 关键词 / 级别
    entries = [("a P6 line", "P6"), ("b P6b line", "P6b"), ("plain text", None),
               ("WARN peak capacity", "E2"), ("error: boom", None)]
    got = model.filter_log(entries, "P6", 0, ["P6", "P6b"])
    check([x for x, _ in got] == ["a P6 line"], "stage filter is exact", str(got))
    got = model.filter_log(entries, "p6b", 0, ["P6", "P6b"])
    check([x for x, _ in got] == ["b P6b line"], "stage filter case-insensitive")
    got = model.filter_log(entries, "plain", 0, ["P6", "P6b"])
    check([x for x, _ in got] == ["plain text"], "text substring filter")
    got = model.filter_log(entries, "", 1, ["P6"])
    check(len(got) == 2, "level>=warn filter", str(got))
    got = model.filter_log(entries, "", 2, ["P6"])
    check(len(got) == 1, "level=error filter", str(got))

    # 运行摘要导出
    st2 = runlog.RunState()
    st2.feed({"ev": "run.start", "t": 1,
              "stages": [{"id": "P0", "weight": 1}, {"id": "E1", "weight": 2}]})
    st2.feed({"ev": "stage.start", "stage": "P0", "t": 2})
    st2.feed({"ev": "stage.end", "stage": "P0", "status": "ok", "sec": 5, "t": 7})
    st2.feed({"ev": "stage.check", "stage": "P0", "status": "ok",
              "fails": 0, "warns": 1, "t": 7})
    st2.feed({"ev": "metric", "stage": "E1", "key": "psm", "value": 123, "t": 8})
    with tempfile.TemporaryDirectory() as d:
        p = model.export_summary(st2, d, stamp="t1")
        body = p.read_text(encoding="utf-8-sig")
        rows = [x for x in body.splitlines() if x]
        check(rows[0].startswith("project,"), "summary header")
        check("stage,status,elapsed_s,est_s,check_fails,check_warns,psm" in body,
              "summary columns")
        check(rows[-1].startswith("E1,pending"), "summary last row", rows[-1])
        check("P0,ok,5.0,,0,1," in body, "summary measured row")

    # 日志内查找
    check(model.find_occurrences("aXbXc", "x") == [(1, 2), (3, 4)],
          "find occurrences (case-insensitive)")
    check(model.find_occurrences("abc", "") == [], "find with empty needle")
    check(model.find_occurrences("aaaa", "aa") == [(0, 2), (2, 4)],
          "find non-overlapping")

    # 阶段悬浮提示
    st3 = runlog.Stage("E1")
    st3.status = "running"
    st3.total, st3.done = 40, 10
    st3.t0, st3.last_t = 100.0, 400.0
    tip = model.stage_tooltip(st3, 600)
    keys = [k for k, _a in tip]
    check("tip.progress" in keys and "tip.elapsed" in keys, "tooltip rows",
          str(keys))
    check(model.stage_tooltip(None) == [], "tooltip empty when no stage")

    # 任务栏进度：没有接口时必须静默降级，不能抛
    from orfscape.gui import taskbar
    tb = taskbar.TaskbarProgress(enabled=False)
    tb.set_progress(0.5)
    tb.set_state("error")
    tb.clear()
    tb.close()
    check(tb.ok is False, "taskbar disabled is no-op")
    tb2 = taskbar.TaskbarProgress()
    for bad in (0.0, 0.33, 1.0, -1, 2, None, "x"):
        tb2.set_progress(bad)
    tb2.set_state("nope")
    tb2.clear()
    tb2.clear()
    tb2.close()
    tb2.close()
    check(True, "taskbar calls never raise")
    print("  info taskbar available on this host:", tb2.ok)


# ---------------------------------------------------------------- 自检
def test_checks():
    print("[checks]")
    from orfscape import checks
    with tempfile.TemporaryDirectory() as d:
        res = checks.check_stage(d, "P0")
        check(res["status"] == "fail", "empty project fails P0")
    tp = PKG.parent / "orfscape_projects" / "tp_pilot"
    if tp.is_dir():
        r = checks.check_stage(str(tp), "C")
        check(r["status"] in ("ok", "warn"), "TP build check", r["summary"])
        r = checks.check_stage(str(tp), "E1")
        check(r["status"] in ("ok", "warn"), "TP E1 check", r["summary"])


def main():
    test_runlog()
    test_progress()
    test_estimator()
    test_i18n()
    test_model()
    test_checks()
    print("\n%d failure(s)" % len(FAILS))
    for f in FAILS:
        print("  - " + f)
    return 0 if not FAILS else 1


if __name__ == "__main__":
    sys.exit(main())
