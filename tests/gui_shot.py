# -*- coding: utf-8 -*-
"""GUI 视觉验收：造一份运行中的状态，截图到 out.png。

    python ORFscape/tests/gui_shot.py [--workdir <proj>] [--lang zh] --out shot.png
"""
import argparse
import os
import sys
import time
from pathlib import Path

HERE = Path(os.path.abspath(__file__)).parent
sys.path.insert(0, str(HERE.parent))

from orfscape import runlog                             # noqa: E402
from orfscape.gui.app import App                        # noqa: E402

ORDER = ["P0", "P1", "D0", "C", "E1", "QC", "E2", "P6", "P6b", "P7",
         "P8a", "P8b", "P8c"]
WEIGHT = {"P0": 1.0, "P1": 3.0, "D0": 10.0, "C": 2.0, "E1": 9.0, "QC": 0.4,
          "E2": 30.0, "P6": 1.5, "P6b": 14.0, "P7": 12.0, "P8a": 0.6,
          "P8b": 4.0, "P8c": 2.0}
EST = {"P0": 180, "P1": 60, "D0": 240, "C": 150, "E1": 1600, "QC": 80,
       "E2": 19600, "P6": 90, "P6b": 7200, "P7": 2400, "P8a": 60, "P8b": 600,
       "P8c": 120}


def demo_events(now):
    t0 = now - 5400
    recs = [{"ev": "run.start", "t": t0, "workdir": "", "raw": "D:/raw",
             "species": "Treponema pallidum", "desc": "TP pilot rerun",
             "lang": "en",
             "stages": [{"id": s, "weight": WEIGHT[s], "est_sec": EST[s]}
                        for s in ORDER],
             "estimate": {"total_sec": sum(EST.values())}}]
    t = t0
    for sid in ("P0", "P1", "D0", "C", "E1", "QC"):
        recs.append({"ev": "stage.start", "stage": sid, "t": t})
        t += EST[sid] * 0.8
        recs.append({"ev": "stage.end", "stage": sid, "status": "ok",
                     "sec": EST[sid] * 0.8, "t": t})
        recs.append({"ev": "stage.check", "stage": sid, "status": "ok",
                     "fails": 0, "warns": 0, "t": t})
    recs.append({"ev": "stage.start", "stage": "E2", "t": t})
    recs.append({"ev": "progress", "stage": "E2", "done": 0, "total": 144,
                 "t": t + 1})
    recs.append({"ev": "progress", "stage": "E2", "done": 61, "total": 144,
                 "t": now})
    recs.append({"ev": "metric", "stage": "E2", "key": "psm", "value": 104233,
                 "t": now})
    lines = [
        (None, "[12:04:11] ORFscape run 20261112_120409  project=orfscape_projects/tp_run2"),
        ("P0", "[12:04:12] >>> P0: python -m orfscape scan --workdir ..."),
        ("P0", "[12:04:12] discovered 36 raw file(s)"),
        ("P0", "[12:04:14] species TP -> Treponema pallidum (taxid=160)"),
        ("P0", "[12:05:30] [stage_check P0] P0: 0 fail, 0 warn"),
        ("E2", "[12:34:02] >>> E2: python -m orfscape e2 --workdir ... --dbs aug,6f"),
        ("E2", "[12:34:05] task count: 144 (raw x CV x db)"),
        ("E2", "[12:41:55]   [8/144] TP_TP_TL_F01_-45__aug OK PSM=8123"),
        ("E2", "[12:52:10]   [24/144] TP_TP_TL_F04_-45__aug OK PSM=7344"),
        ("E2", "[13:06:48]   [48/144] TP_TP_LL_F02_-45__aug OK PSM=9520"),
        ("E2", "[13:22:31]   [64/144] TP_TP_TL_F01_-45__6f OK PSM=2211"),
        ("E2", "[13:38:09]   [61/144] WARN pFind peak capacity 82% — consider more threads"),
    ]
    for i, (stg, line) in enumerate(lines):
        recs.append({"ev": "log", "t": t0 + 60 + i * 420, "line": line,
                     "stage": stg})
    return recs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", default=None)
    ap.add_argument("--lang", default="en")
    ap.add_argument("--tab", default="log")
    ap.add_argument("--theme", default="light")
    ap.add_argument("--scale", type=float, default=1.0)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    if args.workdir is None:
        args.workdir = str(HERE.parent.parent / "orfscape_projects" / "tp_pilot")

    app = App(prefill_workdir=args.workdir)
    app.geometry("1360x744")
    app.set_theme(args.theme, save=False)
    if abs(args.scale - 1.0) > 1e-6:
        app._apply_font_scale(args.scale, save=False)
    app.lang = args.lang
    app.cb_lang.current(1 if args.lang == "zh" else 0)
    app.retranslate()
    app.v_workdir.set(args.workdir)
    app.v_species.set("Treponema pallidum")
    app.v_desc.set("TP pilot · rerun 2")
    app.v_raws = ["D:/raw/XP_SAMPLE_FAIMS_480_DDA_75um_300nL_60min_col0805_20cm_withtrap_L1_45.raw",
                  "D:/raw/XP_SAMPLE_FAIMS_480_DDA_75um_300nL_60min_col0805_20cm_withtrap_L2_45.raw",
                  "D:/raw/XP_SAMPLE_FAIMS_480_DDA_75um_300nL_60min_col0805_20cm_withtrap_L3_45.raw"]
    app._render_rawlist()
    app._refresh_estimate()
    app._refresh_db()

    recs = demo_events(time.time())
    app.state = runlog.RunState()
    for r in recs:
        app.state.feed(r)
        if r.get("ev") == "log":
            app._append_log(r["line"])
    app._apply(recs)

    if args.tab == "db":
        app._open_db_dialog()
    if args.tab == "plan":
        _rows = app._plan_rows(args.workdir, app._form())
        from orfscape import estimator as _est_mod
        _est = _est_mod.estimate(args.workdir, app.v_raws,
                                 threads=app.v_threads.get(),
                                 workers=app.v_workers.get(), dbs=app.v_dbs.get())
        app._plan_win, _choice = app._build_plan_window(_rows, _est)
    app.update_idletasks()
    app.update()
    for _ in range(6):
        app.update()
        time.sleep(0.25)

    app.lift()
    app.attributes('-topmost', True)
    app.update()
    time.sleep(0.35)
    from PIL import ImageGrab
    target = app
    if args.tab == "db" and getattr(app, "_db_win", None) is not None:
        target = app._db_win
    if args.tab == "plan" and getattr(app, "_plan_win", None) is not None:
        target = app._plan_win
    x = target.winfo_rootx()
    y = target.winfo_rooty()
    w = target.winfo_width()
    h = target.winfo_height()
    img = ImageGrab.grab(bbox=(x, y, x + w, y + h))
    img.save(args.out)
    print("saved", args.out, img.size)
    app.destroy()
    return 0


if __name__ == "__main__":
    sys.exit(main())
