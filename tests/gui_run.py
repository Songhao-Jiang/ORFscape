# -*- coding: utf-8 -*-
"""GUI「点开始」链路的无显示验证：运行计划 → 取消分支 → 真拉子进程。

    python ORFscape/tests/gui_run.py [--out shot.png]

断言
- `_plan_rows` 能列出 13 个阶段，首阶段未完成、命令非空；
- 计划确认取消时不启动进程、不写事件文件；
- 确认后子进程被拉起并退出，界面读到 run.start / stage.end / run.end；
- 状态收敛为 ok / fail / check-fail，日志有行，stdout 兜底文件已生成。
"""
import argparse
import os
import sys
import tempfile
import time
from pathlib import Path

HERE = Path(os.path.abspath(__file__)).parent
sys.path.insert(0, str(HERE.parent))

from orfscape.gui.app import App                        # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    fails = []

    def check(cond, label, detail=""):
        if not cond:
            fails.append(label)
            print("FAIL", label, detail)

    with tempfile.TemporaryDirectory() as d:
        proj = os.path.join(d, "proj")
        raw = os.path.join(d, "raw")
        os.makedirs(raw)
        # 占位 .raw：先通过界面"有没有 raw"的预检；解析本身会失败（正是要测的路径）
        (Path(raw) / "dummy.raw").write_bytes(b"RIFF0000")

        app = App()
        app.geometry("1360x744")
        app.v_workdir.set(proj)
        app.v_raws = [raw]
        app._render_rawlist()

        # ---- 1) 运行计划
        rows = app._plan_rows(proj, app._form())
        check(len(rows) == 13, "plan has 13 stages", str(len(rows)))
        check(rows[0][0] == "P0" and not rows[0][2], "plan starts at P0, not done")
        check(all(r[3] and r[3] != "(command unavailable)" for r in rows),
              "every plan row has a command")

        # ---- 2) 取消分支：确认框返回 False 时不应启动
        app._confirm_plan = lambda *a, **k: False
        app.v_confirm.set(True)
        app._start(None)
        app.update()
        ev = os.path.join(proj, "00_logs", "progress.jsonl")
        check(app._proc is None, "cancelled plan did not spawn")
        check(not os.path.exists(ev), "cancelled plan wrote no events")

        # ---- 3) 真跑（跳过确认框）
        app.v_confirm.set(False)
        app._start(None)
        check(app._proc is not None, "process spawned")
        t0 = time.time()
        while time.time() - t0 < 90:
            app.update()
            if app._proc is None:
                break
            time.sleep(0.2)
        app.update()

        st = app.state
        print("status:", st.status, "events:", st.n_events)
        print("stages:", [(s.sid, s.status) for s in st.rows()])
        print("log lines:", len(app._log_data))
        stdout_logs = [f for f in os.listdir(os.path.join(proj, "00_logs"))
                       if f.startswith("gui_stdout_")]
        print("stdout file:", stdout_logs)

        check(app._proc is None, "child exited")
        check(st.status in ("ok", "fail", "check-fail"),
              "terminal status reached", st.status)
        if st.status == "ok":
            check(abs(st.overall() - 1.0) < 1e-6, "ok run is 100%")
        else:
            check(st.failed() is not None, "failure names a stage")
        check(st.n_events >= 3, "events parsed", str(st.n_events))
        check(bool(app._log_data), "log lines rendered")
        check(bool(stdout_logs), "stdout mirror file")
        check(str(app.btn_start["state"]) == "normal", "start button re-enabled")

        # 进度显示：标题栏降级 + 任务栏接口（本机可能没有 ITaskbarList3）
        app._set_title_progress(0.42)
        check("42%" in app.title(), "title shows progress", app.title())
        app._set_title_progress(None)
        check(app.title() == "ORFscape", "title restored", app.title())
        print("taskbar interface available:", app._taskbar.ok)
        app._taskbar.set_progress(0.5)
        app._taskbar.set_state("error")
        app._taskbar.clear()

        # 表格增量刷新：重复渲染不应丢行或抛异常
        before = len(app.tree.get_children())
        for _ in range(100):
            app._render_table()
        check(len(app.tree.get_children()) == before, "incremental table keeps rows",
              "%d -> %d" % (before, len(app.tree.get_children())))
        check(before == 13, "table has 13 stages", str(before))

        if args.out:
            from PIL import ImageGrab
            x, y = app.winfo_rootx(), app.winfo_rooty()
            w, h = app.winfo_width(), app.winfo_height()
            ImageGrab.grab(bbox=(x, y, x + w, y + h)).save(args.out)
            print("saved", args.out)
        app.destroy()

    print("OK" if not fails else "FAILED (%d)" % len(fails))
    return 0 if not fails else 1


if __name__ == "__main__":
    sys.exit(main())
