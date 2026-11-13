# -*- coding: utf-8 -*-
"""ORFscape 图形界面（tkinter）。

入口：
    python -m orfscape gui            # 推荐
    python -m orfscape.gui            # 等价

包本身不 import tkinter：只有真正用界面时才通过 ``__getattr__`` 懒加载
``app``，这样无显示环境（CI、纯命令行批处理）不会因为缺 tk 而失败。
"""


def __getattr__(name):
    if name in ("App", "main"):
        from . import app
        return getattr(app, name)
    raise AttributeError("module %r has no attribute %r" % (__name__, name))


__all__ = ["App", "main"]
