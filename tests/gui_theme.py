# -*- coding: utf-8 -*-
"""主题切换的自动化验收（需要显示环境）。

    python ORFscape/tests/gui_theme.py

断言：
1. 切到 dark 后，界面上不再残留「只在浅色里出现」的颜色值；
2. 切回 light 后，不再残留「只在深色里出现」的颜色值；
3. 主题写进偏好文件（测完还原，不动用户设置）；
4. 两种主题下看板都能重画、日志都能重排，不抛异常。
"""
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

from orfscape.gui import model, theme                 # noqa: E402
from orfscape.gui.app import App                      # noqa: E402

OPTS = ("background", "foreground", "activebackground", "selectbackground",
        "selectforeground", "highlightbackground", "insertbackground",
        "selectcolor", "buttonbackground", "readonlybackground")


def colors_of(w):
    out = {}
    for opt in OPTS:
        try:
            val = w.cget(opt)
        except Exception:
            continue
        if val:
            out[opt] = str(val).strip().lower()
    return out


def leaks(app, values):
    bad = []
    for w in [app] + list(app._iter_widgets()):
        for opt, val in colors_of(w).items():
            if val in values:
                bad.append((w.winfo_class(), opt, val))
    return bad


def main():
    settings = Path(os.path.realpath(os.path.expanduser("~"))) / ".orfscape_gui.json"
    backup = settings.read_bytes() if settings.is_file() else None
    fails = []
    try:
        light_only = {v.lower() for v in theme.LIGHT.values()} - \
                     {v.lower() for v in theme.DARK.values()}
        dark_only = {v.lower() for v in theme.DARK.values()} - \
                    {v.lower() for v in theme.LIGHT.values()}

        app = App()
        app.geometry("1360x744")
        app.set_theme("light")
        app.update()

        app.set_theme("dark")
        app.update()
        bad = leaks(app, light_only)
        if bad:
            fails.append("light colors left in dark mode: %r" % (bad[:6],))
        if model.load_settings().get("theme") != "dark":
            fails.append("theme not persisted")

        app.set_theme("light")
        app.update()
        bad = leaks(app, dark_only)
        if bad:
            fails.append("dark colors left in light mode: %r" % (bad[:6],))
        if model.load_settings().get("theme") != "light":
            fails.append("theme not restored in settings")

        # 主题切换后看板/表格/日志都应可重画
        app._render_board()
        app._render_table()
        app._append_log("theme switch smoke line")
        app.update()
        if not app._log_data:
            fails.append("log lost after theme switch")
        app.destroy()
    finally:
        if backup is not None:
            settings.write_bytes(backup)
        elif settings.is_file():
            settings.unlink()

    if fails:
        for f in fails:
            print("FAIL", f)
        return 1
    print("OK — theme switch is clean (no leaked colors, persisted)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
