# -*- coding: utf-8 -*-
"""GUI「接管已有日志」的无显示验证。

    python ORFscape/tests/gui_attach.py --workdir <project>

流程：造一个项目 + 一份 progress.jsonl（或直接用真跑出来的），让界面
``_attach`` 进去，断言阶段状态、ETA、日志行都到位，然后截图。
"""
import argparse
import os
import sys
import time
from pathlib import Path

HERE = Path(os.path.abspath(__file__)).parent
sys.path.insert(0, str(HERE.parent))

from orfscape.gui.app import App                        # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    app = App(prefill_workdir=args.workdir)
    app.geometry("1360x744")
    ok = app._attach(args.workdir, quiet=True)
    app.update_idletasks()
    app.update()
    st = app.state
    print("attached:", ok)
    print("events:", st.n_events, "status:", st.status)
    print("stages:", [(s.sid, s.status) for s in st.rows()])
    print("overall: %.3f" % st.overall())
    print("eta:", st.eta_sec())
    print("log lines:", len(app._log_data))
    failed = [s.sid for s in st.rows() if s.status in ("fail", "check-fail")]
    if ok:
        assert st.n_events > 0, "no events read"
        assert app._log_data, "no log lines rendered"
        print("failed stages:", failed)
    if args.out:
        for _ in range(4):
            app.update()
            time.sleep(0.2)
        from PIL import ImageGrab
        x, y = app.winfo_rootx(), app.winfo_rooty()
        w, h = app.winfo_width(), app.winfo_height()
        ImageGrab.grab(bbox=(x, y, x + w, y + h)).save(args.out)
        print("saved", args.out)
    app.destroy()
    return 0


if __name__ == "__main__":
    sys.exit(main())
