# -*- coding: utf-8 -*-
"""界面视觉回归：渲染固定状态，与基线图比像素差。

    python ORFscape/tests/gui_visual.py                 # 比对（无基线则创建）
    python ORFscape/tests/gui_visual.py --update        # 重建基线
    python ORFscape/tests/gui_visual.py --threshold 6

只比较「下半部分」（resources 条以下），避开实时 CPU/内存数字；其余区域在固定
状态+固定窗口尺寸下应当是确定的。基线放在 `tests/baseline/gui_state.png`，
换机器/换字体/改主题后需要 `--update` 重录。
"""
import argparse
import os
import sys
import tempfile
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))

from gui_shot import demo_events                     # noqa: E402
from orfscape import runlog                            # noqa: E402
from orfscape.gui.app import App                       # noqa: E402

BASELINE = HERE / "baseline" / "gui_state.png"
CROP_TOP = 150          # 跳过 header + 资源条
CROP_BOTTOM = 70        # 跳过状态栏


def shoot(path, workdir):
    app = App(prefill_workdir=workdir)
    app.geometry("1360x744")
    app.set_theme("light", save=False)
    app.v_workdir.set(workdir)
    app.v_species.set("Treponema pallidum")
    app.v_desc.set("TP pilot · rerun 2")
    app.v_raws = ["D:/raw/XP_SAMPLE_FAIMS_480_DDA_75um_300nL_60min_col0805_20cm_withtrap_L1_45.raw"]
    app._render_rawlist()
    app._refresh_db()
    recs = demo_events(time.time())
    app.state = runlog.RunState()
    for r in recs:
        app.state.feed(r)
        if r.get("ev") == "log":
            app._append_log(r["line"])
    app._apply(recs)
    app.update_idletasks()
    for _ in range(5):
        app.update()
        time.sleep(0.15)
    app.lift()
    app.attributes('-topmost', True)
    app.update()
    time.sleep(0.35)
    from PIL import ImageGrab
    x, y = app.winfo_rootx(), app.winfo_rooty()
    w, h = app.winfo_width(), app.winfo_height()
    ImageGrab.grab(bbox=(x, y, x + w, y + h)).save(path)
    app.destroy()


def compare(a, b):
    from PIL import Image, ImageChops
    ia, ib = Image.open(a).convert("RGB"), Image.open(b).convert("RGB")
    if ia.size != ib.size:
        return None, "size %s != %s" % (ia.size, ib.size)
    box = (0, CROP_TOP, min(ia.size[0], ib.size[0]),
           min(ia.size[1], ib.size[1]) - CROP_BOTTOM)
    ca, cb = ia.crop(box), ib.crop(box)
    diff = ImageChops.difference(ca, cb)
    px = list(diff.getdata())
    total = sum(sum(p) for p in px) / float(len(px) * 3)
    return total, box


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--update", action="store_true")
    ap.add_argument("--threshold", type=float, default=5.0)
    ap.add_argument("--workdir", default=None)
    args = ap.parse_args()
    workdir = args.workdir or str(HERE.parent.parent / "orfscape_projects" / "tp_pilot")

    with tempfile.TemporaryDirectory() as d:
        shot = os.path.join(d, "shot.png")
        shoot(shot, workdir)
        if args.update or not BASELINE.is_file():
            os.makedirs(BASELINE.parent, exist_ok=True)
            from PIL import Image
            Image.open(shot).save(BASELINE)
            print("baseline written:", BASELINE)
            return 0
        score, box = compare(shot, str(BASELINE))
        if score is None:
            print("FAIL", box)
            return 1
        print("mean abs diff = %.2f (threshold %.2f, box %s)"
              % (score, args.threshold, box))
        if score > args.threshold:
            print("FAIL visual regression — rerun with --update if the change is intended")
            return 1
        print("OK")
        return 0


if __name__ == "__main__":
    sys.exit(main())
