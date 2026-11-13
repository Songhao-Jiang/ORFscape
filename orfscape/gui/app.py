# -*- coding: utf-8 -*-
"""ORFscape 的图形界面（tkinter，单窗口）。

界面语言默认英文，右上角可切中文；所有可见文案走 :mod:`orfscape.gui.i18n`。

设计
----
- **界面只读日志**：流程由 ``python -m orfscape run`` 子进程跑，它把结构化事件
  写进 ``<项目>/00_logs/progress.jsonl``；界面轮询这个文件更新进度、ETA、
  阶段状态与日志，因此同一个文件也能被 agent/skill 的跑法复用。
- **估算与实测分离**：开跑前用 :mod:`orfscape.estimator` 给量级（标注 estimated
  / measured），跑起来后 ETA 由实测速率接管。
- **资源**：:mod:`orfscape.sysmon` 直接问内核（Windows ctypes / Linux /proc），
  不依赖 psutil；取不到就显示 n/a。
"""
import datetime
import os
import queue
import re
import subprocess
import sys
import threading
import types
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

from .. import estimator, runlog, sysmon
from ..estimator import fmt_bytes, fmt_sec
from . import i18n, model, taskbar, theme
from .theme import C


def scaled_font(size, weight="normal", scale=1.0):
    """按缩放系数给出字体元组（canvas 文本用）。"""
    return ("Segoe UI", max(6, int(round(size * scale))), weight)


def st_style(status):
    """阶段状态 → (边框色, 底色, 圆点色)；跟随主题实时取色。"""
    if status == "running":
        return C["accent"], C["chip_run"], C["accent"]
    if status == "ok":
        return C["ok_border"], C["chip_ok"], C["ok"]
    if status in ("fail", "check-fail"):
        return C["bad_border"], C["chip_bad"], C["fail"]
    if status == "stopped":
        return C["stop_border"], C["stop_bg"], C["warn"]
    return C["border"], C["chip"], C["idle"]

FONT = ("Segoe UI", 9)
FONT_B = ("Segoe UI", 9, "bold")
FONT_H2 = ("Segoe UI", 10, "bold")
FONT_BIG = ("Segoe UI", 14, "bold")
FONT_MONO = ("Consolas", 9)


def _rr(canvas, x0, y0, x1, y1, r, **kw):
    """圆角矩形（多边形近似）。"""
    r = max(0, min(r, (x1 - x0) / 2, (y1 - y0) / 2))
    pts = [x0 + r, y0, x1 - r, y0, x1, y0, x1, y0 + r, x1, y1 - r, x1, y1,
           x1 - r, y1, x0 + r, y1, x0, y1, x0, y1 - r, x0, y0 + r, x0, y0]
    return canvas.create_polygon(pts, smooth=True, **kw)


class App(tk.Tk):
    # ------------------------------------------------------------ 初始化
    def __init__(self, prefill_workdir=None):
        tk.Tk.__init__(self)
        self._dpi()
        self.cfg = model.load_settings()
        self.lang = self.cfg.get("lang") or i18n.LANG_EN
        self.state = runlog.RunState()
        self._ev_off = 0
        self._events_path = None
        self._proc = None
        self._stdout_fh = None
        self._run_t0 = None
        self._stopping = False
        self._bindings = []
        self._hooks = []
        self._db_rows = []
        self._est = None
        self._sampler = sysmon.Sampler()
        self._last_sample = 0.0
        self._sample = {}
        self._log_lines = 0
        self._ui_q = queue.Queue()
        self._title_pct = None
        self._font_scale = float(self.cfg.get("font_scale") or 1.0)
        self._font_base = {}
        self._taskbar = taskbar.TaskbarProgress()

        self.title("ORFscape")
        theme.set_theme(self.cfg.get("theme") or "light")
        sw, sh = self.winfo_screenwidth(), self.winfo_screenheight()
        self.geometry("%dx%d" % (min(1400, sw - 60), min(900, sh - 120)))
        self.minsize(1120, 700)
        self.configure(bg=C["bg"])
        self._style()
        self._vars()
        self._build()
        self._taskbar.bind(self)
        if abs(self._font_scale - 1.0) > 1e-6:
            self._apply_font_scale(self._font_scale, save=False)
        self._restore(prefill_workdir)
        self.retranslate()
        self.after(500, self._poll)
        self.after(200, self._sample_resources)
        self.after(700, self._refresh_estimate)
        self.protocol("WM_DELETE_WINDOW", self._on_close)
        geo = str(self.cfg.get("geometry") or "")
        if re.match(r"^\d+x\d+([+-]-?\d+){1,2}$", geo):
            try:
                self.geometry(geo)
            except tk.TclError:
                pass
        self._bind_keys()

    def _bind_keys(self):
        """快捷键：Ctrl+Enter 开始 / Esc 停止 / Ctrl+L 过滤框 / F5 重算。"""
        self.bind_all("<Control-Return>",
                      lambda e: (self._start(None), "break")[1])
        self.bind_all("<Escape>", self._esc)
        self.bind_all("<Control-l>",
                      lambda e: (self.ent_filter.focus_set(), "break")[1])
        self.bind_all("<F5>", lambda e: (self._refresh_estimate(),
                                         self._refresh_db(), "break")[2])
        self.bind_all("<Control-f>", lambda e: self._toggle_find())
        self.bind_all("<Control-o>", lambda e: (self._add_files(), "break")[1])
        self.bind_all("<Control-d>", lambda e: (self._open_db_dialog(), "break")[1])
        self.bind_all("<Control-Shift-C>", lambda e: (self._copy_all(), "break")[1])

    def _esc(self, _ev=None):
        try:
            cls = self.focus_get().winfo_class()
        except (AttributeError, tk.TclError):
            cls = ""
        if cls in ("Entry", "TEntry", "Text", "Listbox", "TCombobox",
                   "Spinbox", "TSpinbox"):
            return None
        if self._proc is not None:
            self._stop()
        return "break"

    def _dpi(self):
        if os.name == "nt":
            try:
                import ctypes
                ctypes.windll.shcore.SetProcessDpiAwareness(1)
            except Exception:
                pass

    def _style(self):
        st = ttk.Style(self)
        try:
            st.theme_use("clam")
        except tk.TclError:
            pass
        st.configure(".", font=FONT, background=C["bg"], foreground=C["text"])
        st.configure("Card.TFrame", background=C["card"])
        st.configure("Bar.Horizontal.TProgressbar", troughcolor=C["track"],
                     background=C["accent"], bordercolor=C["track"],
                     lightcolor=C["accent"], darkcolor=C["accent"], thickness=8)
        st.configure("Slim.Horizontal.TProgressbar", troughcolor=C["track"],
                     background=C["accent"], bordercolor=C["track"],
                     lightcolor=C["accent"], darkcolor=C["accent"], thickness=5)
        st.configure("Treeview", background=C["card"], fieldbackground=C["card"],
                     foreground=C["text"], rowheight=21, borderwidth=0)
        st.configure("Treeview.Heading", font=FONT_B, background=C["bg"],
                     foreground=C["muted"], relief="flat")
        st.map("Treeview", background=[("selected", C["sel"])],
               foreground=[("selected", C["text"])])
        st.configure("TEntry", fieldbackground=C["field"], foreground=C["text"],
                     bordercolor=C["border"], insertcolor=C["text"],
                     lightcolor=C["border"], darkcolor=C["border"])
        st.configure("TCombobox", fieldbackground=C["field"], foreground=C["text"],
                     arrowcolor=C["muted"], bordercolor=C["border"],
                     background=C["chip"], selectbackground=C["sel"],
                     selectforeground=C["text"], lightcolor=C["border"],
                     darkcolor=C["border"])
        # clam 自带 readonly 状态映射，会盖掉上面的 configure，必须显式覆盖
        st.map("TCombobox",
               fieldbackground=[("readonly", C["field"]), ("disabled", C["bg"])],
               foreground=[("readonly", C["text"]), ("disabled", C["idle"])],
               selectbackground=[("readonly", C["sel"])],
               selectforeground=[("readonly", C["text"])],
               arrowcolor=[("disabled", C["idle"])],
               background=[("active", C["chip"]), ("pressed", C["chip"]),
                           ("readonly", C["chip"])])
        # 下拉列表（Tk listbox）不受 ttk 样式影响，走 option 数据库
        self.option_add("*TCombobox*Listbox.background", C["field"])
        self.option_add("*TCombobox*Listbox.foreground", C["text"])
        self.option_add("*TCombobox*Listbox.selectBackground", C["sel"])
        self.option_add("*TCombobox*Listbox.selectForeground", C["text"])
        for sb in ("Vertical.TScrollbar", "Horizontal.TScrollbar", "TScrollbar"):
            st.configure(sb, background=C["border"], troughcolor=C["bg"],
                         bordercolor=C["bg"], arrowcolor=C["muted"],
                         lightcolor=C["border"], darkcolor=C["border"],
                         relief="flat")
        st.configure("TButton", background=C["chip"], foreground=C["text"],
                     bordercolor=C["border"], focusthickness=1,
                     focuscolor=C["accent"], relief="flat", padding=(8, 3))
        st.map("TButton", background=[("active", C["sel"]),
                                      ("disabled", C["chip"])],
               foreground=[("disabled", C["idle"])])
        st.configure("TNotebook", background=C["card"], bordercolor=C["border"])
        st.configure("TNotebook.Tab", background=C["chip"], foreground=C["muted"],
                     padding=(10, 3))
        st.map("TNotebook.Tab", background=[("selected", C["card"])],
               foreground=[("selected", C["text"])])

    # ------------------------------------------------------------ 字号缩放
    def _fnt(self, size, weight="normal"):
        return scaled_font(size, weight, getattr(self, "_font_scale", 1.0))

    def _apply_font_scale(self, scale, save=True):
        """按系数重设所有控件的字号与表格行高（基准值只在首次记录）。"""
        import tkinter.font as tkfont
        try:
            scale = float(scale)
        except (TypeError, ValueError):
            return
        scale = min(1.6, max(0.8, scale))
        base_cache = self._font_base
        for w in [self] + list(self._iter_widgets()):
            key = str(w)
            try:
                cur = w.cget("font")
            except tk.TclError:
                continue
            if not cur:
                continue
            spec = base_cache.get(key)
            if spec is None:
                try:
                    a = tkfont.Font(font=cur).actual()
                except Exception:
                    continue
                spec = (a.get("family") or "Segoe UI", int(a.get("size") or 9),
                        a.get("weight") or "normal", a.get("slant") or "roman")
                base_cache[key] = spec
            size = max(6, int(round(spec[1] * scale)))
            try:
                w.configure(font=(spec[0], size,
                                  "bold" if spec[2] == "bold" else "normal",
                                  "italic" if spec[3] == "italic" else "roman"))
            except tk.TclError:
                pass
        self._font_scale = scale
        if hasattr(self, "_hdr"):
            self._hdr.configure(height=int(round(56 * scale)))
        if hasattr(self, "board"):
            # 看板不随字号无限长高，否则把日志挤没
            self.board.configure(height=min(150, int(round(112 * scale))))
        if hasattr(self, "tree"):
            # 字号越大，明细表行数越少（看板才是主视图）
            self.tree.configure(
                height=max(3, int(round(5 - (scale - 1.0) * 7))))
        st = ttk.Style(self)
        st.configure("Treeview", rowheight=max(16, int(round(21 * scale))))
        st.configure("Treeview.Heading",
                     font=(FONT_B[0], max(7, int(round(FONT_B[1] * scale))),
                           "bold"))
        if hasattr(self, "tree"):
            widths = {"#0": 190, "status": 90, "progress": 110, "elapsed": 80,
                      "eta": 80, "check": 110}
            for col, wd in widths.items():
                self.tree.column(col, width=int(wd * scale))
        self._render_board()
        self._render_table()
        if save:
            self.cfg["font_scale"] = scale
            self._save_cfg()

    # ------------------------------------------------------------ 主题
    def _toggle_theme(self):
        self.set_theme("light" if theme.get_theme() == "dark" else "dark")

    def _iter_widgets(self, root=None):
        for child in (root or self).winfo_children():
            yield child
            for sub in self._iter_widgets(child):
                yield sub

    def _recolor_tree(self, mapping):
        """按"旧色→新色"映射重涂 tk 控件的颜色选项（ttk 交给 _style）。"""
        opts = ("background", "foreground", "activebackground",
                "activeforeground", "selectbackground", "selectforeground",
                "highlightbackground", "highlightcolor", "insertbackground",
                "selectcolor", "troughcolor", "buttonbackground",
                "readonlybackground")
        for w in [self] + list(self._iter_widgets()):
            for opt in opts:
                try:
                    cur = w.cget(opt)
                except tk.TclError:
                    continue
                if not cur:
                    continue
                new = theme.map_color(cur, mapping)
                if new == cur:
                    continue
                try:
                    w.configure(**{opt: new})
                except tk.TclError:
                    pass

    def _apply_roles(self):
        """同色 token 的角色差异（输入框用 field、日志区用 log_bg…）。"""
        for w in [self] + list(self._iter_widgets()):
            cls = w.winfo_class()
            try:
                if cls in ("Entry", "Spinbox"):
                    w.configure(bg=C["field"], fg=C["text"],
                                insertbackground=C["text"])
                elif cls == "Text":
                    if w is self.txt:
                        w.configure(bg=C["log_bg"], fg=C["log_fg"])
                    else:
                        w.configure(bg=C["field"], fg=C["text"])
                elif cls == "Listbox":
                    w.configure(bg=C["field"], fg=C["text"],
                                selectbackground=C["sel"],
                                selectforeground=C["text"])
            except tk.TclError:
                pass
        self.txt.tag_configure("warn", foreground=C["warn"])
        self.txt.tag_configure("err", foreground=C["fail"])
        self.txt.tag_configure("ok", foreground=C["ok"])
        self.txt.tag_configure("dim", foreground=C["hint"])
        self.txt.tag_configure("hit", background=C["hit"])
        self.txt.tag_configure("find", background=C["find_bg"])
        self.txt.tag_configure("find_cur", background=C["accent"],
                               foreground=C["field"])
        self.tree.tag_configure("running", foreground=C["accent"],
                                background=C["chip_run"])
        self.tree.tag_configure("fail", foreground=C["fail"],
                                background=C["chip_bad"])
        self.tree.tag_configure("check-fail", foreground=C["fail"],
                                background=C["chip_bad"])
        self.tree.tag_configure("ok", foreground=C["text"])
        self.tree.tag_configure("pending", foreground=C["muted"])
        self.tree.tag_configure("odd", background=C["zebra"])
        self._render_theme_label()

    def _render_theme_label(self):
        if not hasattr(self, "btn_theme"):
            return
        key = "theme.to_light" if theme.get_theme() == "dark" else "theme.to_dark"
        self.btn_theme.configure(text=self.t(key))

    def set_theme(self, name, save=True):
        """切换明暗主题：换调色板 → 重涂 → 重建样式 → 重画。

        save=False 用于截图/测试，不写用户偏好。
        """
        if name == theme.get_theme():
            self._render_theme_label()
            return
        mapping = theme.color_map()
        theme.set_theme(name)
        self.configure(bg=C["bg"])
        self._recolor_tree(mapping)
        self._apply_roles()
        self._style()
        self._render_board()
        self._render_table()
        self._render_db()
        self._render_log(force=True)
        if save:
            self.cfg["theme"] = theme.get_theme()
            self._save_cfg()

    def _vars(self):
        self.v_raws = []
        self.v_species = tk.StringVar(value=self.cfg.get("species", ""))
        self.v_desc = tk.StringVar(value=self.cfg.get("desc", ""))
        self.v_workdir = tk.StringVar(value=self.cfg.get("workdir", ""))
        self.v_threads = tk.IntVar(value=int(self.cfg.get("threads") or 16))
        self.v_workers = tk.IntVar(value=int(self.cfg.get("workers") or 4))
        self.v_dbs = tk.StringVar(value=self.cfg.get("dbs") or "aug,6f")
        self.v_skip_p6b = tk.BooleanVar(value=bool(self.cfg.get("skip_p6b")))
        self.v_skip_qc = tk.BooleanVar(value=bool(self.cfg.get("skip_qc")))
        self.v_force = tk.BooleanVar(value=bool(self.cfg.get("force")))
        self.v_confirm = tk.BooleanVar(value=bool(self.cfg.get("confirm", True)))
        self.v_lang = tk.StringVar(value=self.lang)
        self.v_filter = tk.StringVar(value=self.cfg.get("filter") or "")
        self.v_autoscroll = tk.BooleanVar(value=True)
        self.v_status = tk.StringVar(value="")
        self.v_current = tk.StringVar(value="")
        self.v_elapsed = tk.StringVar(value="")
        self.v_eta = tk.StringVar(value="")
        self.v_overall = tk.StringVar(value="0%")
        self.v_rawsum = tk.StringVar(value="")
        self.v_db_total = tk.StringVar(value="")
        self.v_meta = tk.StringVar(value="")
        self._level_idx = int(self.cfg.get("log_level") or 0)  # 0=all 1=warn+ 2=error

    # ------------------------------------------------------------ 绑定/翻译
    def t(self, key, **kw):
        return i18n.tr(key, self.lang, **kw)

    def _bind_text(self, widget, option, key):
        self._bindings.append((widget, option, key))
        try:
            widget.configure(**{option: self.t(key)})
        except tk.TclError:
            pass

    def _hook(self, fn):
        self._hooks.append(fn)

    def retranslate(self):
        for widget, option, key in self._bindings:
            try:
                widget.configure(**{option: self.t(key)})
            except tk.TclError:
                pass
        for pane, key in getattr(self, "_nb_tabs", []):
            try:
                self.nb.tab(pane, text=self.t(key))
            except tk.TclError:
                pass
        for fn in self._hooks:
            try:
                fn()
            except Exception:
                pass
        if hasattr(self, "btn_adv"):
            self._render_adv_label()
        self._render_level_values()
        self._render_theme_label()
        self._render_log_menu_labels()
        self._render_board()
        self._render_table_headings()
        self._render_db()
        self._refresh_rawsum()

    # ------------------------------------------------------------ 布局
    def _build(self):
        root = tk.Frame(self, bg=C["bg"])
        root.pack(fill="both", expand=True)
        self._build_header(root)
        self._build_resource_strip(root)
        self._build_statusbar(root)          # 先占底部，body 再吃剩余空间

        body = tk.Frame(root, bg=C["bg"])
        body.pack(fill="both", expand=True, padx=10, pady=(0, 6))
        body.grid_columnconfigure(0, weight=0, minsize=430)
        body.grid_columnconfigure(1, weight=1)
        body.grid_rowconfigure(0, weight=1)

        left = tk.Frame(body, bg=C["bg"], width=430)
        left.grid(row=0, column=0, sticky="nsew", padx=(0, 8))
        left.grid_propagate(False)
        self._left_scroll = tk.Canvas(left, bg=C["bg"], highlightthickness=0)
        left_sb = ttk.Scrollbar(left, orient="vertical",
                                command=self._left_scroll.yview)
        self._left_scroll.configure(yscrollcommand=left_sb.set)
        left_sb.pack(side="right", fill="y")
        self._left_scroll.pack(side="left", fill="both", expand=True)
        left_inner = tk.Frame(self._left_scroll, bg=C["bg"])
        self._left_win = self._left_scroll.create_window(
            (0, 0), window=left_inner, anchor="nw")
        left_inner.bind(
            "<Configure>",
            lambda e: self._left_scroll.configure(
                scrollregion=self._left_scroll.bbox("all")))
        self._left_scroll.bind(
            "<Configure>",
            lambda e: self._left_scroll.itemconfigure(self._left_win,
                                                      width=e.width))
        left_inner.bind("<Enter>", lambda e: self._bind_wheel(True))
        left_inner.bind("<Leave>", lambda e: self._bind_wheel(False))
        self._build_input(left_inner)
        self._build_estimate(left_inner)

        right = tk.Frame(body, bg=C["bg"])
        right.grid(row=0, column=1, sticky="nsew")
        right.grid_rowconfigure(0, weight=0)
        right.grid_rowconfigure(1, weight=1, minsize=140)
        right.grid_columnconfigure(0, weight=1)
        top_pane = tk.Frame(right, bg=C["bg"])
        top_pane.grid(row=0, column=0, sticky="nsew")
        self._build_pipeline(top_pane)

        log_pane = tk.Frame(right, bg=C["bg"])
        log_pane.grid(row=1, column=0, sticky="nsew", pady=(8, 0))
        self._build_log(log_pane)

    def _bind_wheel(self, on):
        """左栏滚轮：只在指针位于左栏时接管，避免抢日志/列表的滚轮。"""
        if on:
            self.bind_all("<MouseWheel>", self._on_wheel)
        else:
            self.unbind_all("<MouseWheel>")

    def _on_wheel(self, ev):
        try:
            self._left_scroll.yview_scroll(int(-1 * (ev.delta / 120.0)), "units")
        except tk.TclError:
            pass

    def _build_header(self, root):
        h = tk.Frame(root, bg=C["card"], height=56)
        self._hdr = h
        h.pack(fill="x", padx=10, pady=(5, 5))
        h.pack_propagate(False)
        left = tk.Frame(h, bg=C["card"])
        left.pack(side="left", padx=14, pady=5)
        ttl = tk.Label(left, text="", bg=C["card"], fg=C["text"], font=FONT_BIG)
        ttl.pack(anchor="w")
        self._bind_text(ttl, "text", "app.title")
        sub = tk.Label(left, text="", bg=C["card"], fg=C["muted"], font=FONT)
        sub.pack(anchor="w")
        self._bind_text(sub, "text", "app.subtitle")

        right = tk.Frame(h, bg=C["card"])
        right.pack(side="right", padx=14)
        lab = tk.Label(right, text="", bg=C["card"], fg=C["muted"], font=FONT)
        lab.pack(side="left", padx=(0, 6))
        self._bind_text(lab, "text", "app.lang")
        self.btn_theme = ttk.Button(right, text="", width=11,
                                    command=self._toggle_theme)
        self.btn_theme.pack(side="left", padx=(0, 12))
        db_btn = ttk.Button(right, text="", command=self._open_db_dialog)
        db_btn.pack(side="left", padx=(0, 16))
        self._bind_text(db_btn, "text", "db.open")

        self.cb_lang = ttk.Combobox(right, width=9, state="readonly",
                                    values=[n for _c, n in i18n.LANGS])
        self.cb_lang.current(0 if self.lang == i18n.LANG_EN else 1)
        self.cb_lang.pack(side="left")
        self.cb_lang.bind("<<ComboboxSelected>>", self._on_lang)

    def _card(self, parent, title_key, pady=(0, 6), fill="both", expand=True):
        outer = tk.Frame(parent, bg=C["border"])
        outer.pack(fill=fill, expand=expand, pady=pady)
        inner = tk.Frame(outer, bg=C["card"])
        inner.pack(fill="both", expand=True, padx=1, pady=1)
        head = tk.Frame(inner, bg=C["card"])
        head.pack(fill="x", padx=12, pady=(9, 2))
        lab = tk.Label(head, text="", bg=C["card"], fg=C["text"], font=FONT_H2)
        lab.pack(side="left")
        self._bind_text(lab, "text", title_key)
        body = tk.Frame(inner, bg=C["card"])
        body.pack(fill="both", expand=True, padx=12, pady=(2, 10))
        return outer, body, head

    # ------------------------------------------------------------ 输入卡
    def _build_input(self, parent):
        _o, body, head = self._card(parent, "sec.input", pady=(0, 6),
                                    fill="x", expand=False)
        bar = tk.Frame(body, bg=C["card"])
        bar.pack(fill="x", pady=(0, 6))
        for key, fn, w in (("inp.add_files", self._add_files, 13),
                           ("inp.add_folder", self._add_folder, 13)):
            b = ttk.Button(bar, text="", command=fn, width=w)
            b.pack(side="left", padx=(0, 4))
            self._bind_text(b, "text", key)
        for key, fn, w in (("inp.attach", self._load_existing, 13),
                           ("inp.clear", self._clear_raws, 8)):
            b = ttk.Button(bar, text="", command=fn, width=w)
            b.pack(side="left", padx=(0, 4))
            self._bind_text(b, "text", key)

        wrap = tk.Frame(body, bg=C["border"])
        wrap.pack(fill="x", pady=(2, 3))
        self.lb_raw = tk.Listbox(wrap, height=3, activestyle="none", bd=0,
                                 font=FONT, bg=C["field"], fg=C["text"],
                                 selectbackground=C["sel"],
                                 selectforeground=C["text"], highlightthickness=0,
                                 exportselection=False)
        self.lb_raw.pack(side="left", fill="both", expand=True, padx=1, pady=1)
        sb = ttk.Scrollbar(wrap, orient="vertical", command=self.lb_raw.yview)
        sb.pack(side="right", fill="y")
        self.lb_raw.configure(yscrollcommand=sb.set)
        self.lb_raw.bind("<Delete>", lambda e: self._remove_selected())

        row = tk.Frame(body, bg=C["card"])
        row.pack(fill="x", pady=(0, 4))
        self.lb_sum = tk.Label(row, textvar=self.v_rawsum, bg=C["card"],
                               fg=C["muted"], font=FONT, anchor="w")
        self.lb_sum.pack(side="left")
        b = ttk.Button(row, text="", width=7, command=self._remove_selected)
        b.pack(side="right")
        self._bind_text(b, "text", "inp.remove")

        two = tk.Frame(body, bg=C["card"])
        two.pack(fill="x", pady=(2, 0))
        two.grid_columnconfigure(0, weight=3)
        two.grid_columnconfigure(1, weight=2)
        self._field(two, "inp.species", self.v_species, None, col=0)
        self._field(two, "inp.desc", self.v_desc, None, col=1)
        hint = tk.Label(body, text="", bg=C["card"], fg=C["hint"],
                        font=("Segoe UI", 8), anchor="w")
        hint.pack(fill="x")
        self._bind_text(hint, "text", "inp.species_hint")

        wd = tk.Frame(body, bg=C["card"])
        wd.pack(fill="x", pady=(4, 2))
        lab = tk.Label(wd, text="", bg=C["card"], fg=C["muted"], font=FONT)
        lab.pack(anchor="w")
        self._bind_text(lab, "text", "inp.workdir")
        row = tk.Frame(wd, bg=C["card"])
        row.pack(fill="x")
        ttk.Entry(row, textvariable=self.v_workdir).pack(
            side="left", fill="x", expand=True)
        b = ttk.Button(row, text="", width=9, command=self._pick_workdir)
        b.pack(side="left", padx=(4, 0))
        self._bind_text(b, "text", "inp.browse")

        # 高级项默认收起：小屏（864 px）也能一屏看到估算
        self._adv_open = tk.BooleanVar(value=bool(self.cfg.get("advanced")))
        self.btn_adv = tk.Label(body, text="", bg=C["card"], fg=C["accent"],
                                font=("Segoe UI", 8), anchor="w", cursor="hand2")
        self.btn_adv.pack(fill="x", pady=(6, 0))
        self.btn_adv.bind("<Button-1>", lambda e: self._toggle_adv())

        self.adv = tk.Frame(body, bg=C["card"])
        self._adv_built = False
        if self._adv_open.get():
            self._build_adv()
            self.adv.pack(fill="x", pady=(4, 0))

    def _toggle_adv(self):
        if not self._adv_built:
            self._build_adv()
        if self._adv_open.get():
            self.adv.pack_forget()
            self._adv_open.set(False)
        else:
            self.adv.pack(fill="x", pady=(4, 0))
            self._adv_open.set(True)
        self._render_adv_label()

    def _render_adv_label(self):
        arrow = "▾" if self._adv_open.get() else "▸"
        self.btn_adv.configure(text="%s %s" % (arrow, self.t("inp.advanced")))

    def _build_adv(self):
        row2 = tk.Frame(self.adv, bg=C["card"])
        row2.pack(fill="x", pady=(0, 4))
        for key, var, lo, hi in (("inp.threads", self.v_threads, 1, 128),
                                 ("inp.workers", self.v_workers, 1, 32)):
            lab = tk.Label(row2, text="", bg=C["card"], fg=C["muted"], font=FONT)
            lab.pack(side="left")
            self._bind_text(lab, "text", key)
            tk.Spinbox(row2, from_=lo, to=hi, width=5, textvariable=var,
                       font=FONT, bd=0, highlightthickness=1, bg=C["field"],
                       fg=C["text"], buttonbackground=C["chip"],
                       insertbackground=C["text"],
                       highlightbackground=C["border"]).pack(side="left", padx=(4, 12))

        lab = tk.Label(row2, text="", bg=C["card"], fg=C["muted"], font=FONT)
        lab.pack(side="left")
        self._bind_text(lab, "text", "inp.dbs")
        ttk.Entry(row2, textvariable=self.v_dbs, width=8).pack(side="left", padx=4)

        lab = tk.Label(row2, text="", bg=C["card"], fg=C["muted"], font=FONT)
        lab.pack(side="left", padx=(8, 0))
        self._bind_text(lab, "text", "inp.font_scale")
        self.cb_font = ttk.Combobox(row2, width=6, state="readonly",
                                    values=("90%", "100%", "115%", "130%"))
        self.cb_font.pack(side="left", padx=4)
        opts = (0.9, 1.0, 1.15, 1.3)
        best = min(range(len(opts)), key=lambda i: abs(opts[i] - self._font_scale))
        self.cb_font.current(best)
        self.cb_font.bind("<<ComboboxSelected>>", self._on_font_scale)

        for key, var in (("inp.confirm", self.v_confirm),
                         ("inp.skip_p6b", self.v_skip_p6b),
                         ("inp.skip_qc", self.v_skip_qc),
                         ("inp.force", self.v_force)):
            cb = tk.Checkbutton(self.adv, text="", variable=var, bg=C["card"],
                                fg=C["text"], font=FONT, anchor="w",
                                activebackground=C["card"],
                                selectcolor=C["field"], bd=0,
                                highlightthickness=0)
            cb.pack(anchor="w")
            self._bind_text(cb, "text", key)
        self._adv_built = True

    def _field(self, parent, label_key, var, hint_key=None, col=None):
        wrap = tk.Frame(parent, bg=C["card"])
        if col is None:
            wrap.pack(fill="x", pady=(2, 2))
        else:
            wrap.grid(row=0, column=col, sticky="ew",
                      padx=(0, 8) if col == 0 else (0, 0))
        lab = tk.Label(wrap, text="", bg=C["card"], fg=C["muted"], font=FONT)
        lab.pack(anchor="w")
        self._bind_text(lab, "text", label_key)
        ttk.Entry(wrap, textvariable=var).pack(fill="x")
        if hint_key:
            hint = tk.Label(wrap, text="", bg=C["card"], fg=C["hint"], font=("Segoe UI", 8))
            hint.pack(anchor="w")
            self._bind_text(hint, "text", hint_key)

    # ------------------------------------------------------------ 估算卡
    def _build_estimate(self, parent):
        _o, body, head = self._card(parent, "est.title", fill="x", expand=False)
        b = ttk.Button(head, text="", width=13, command=self._refresh_estimate)
        b.pack(side="right")
        self._bind_text(b, "text", "est.recalc")

        grid = tk.Frame(body, bg=C["card"])
        grid.pack(fill="x")
        grid.grid_columnconfigure(1, weight=1)
        self.est_labels = {}
        for i, key in enumerate(("est.raw", "est.ms2", "est.genome", "est.tasks",
                                 "est.time", "est.disk", "est.ram")):
            lab = tk.Label(grid, text="", bg=C["card"], fg=C["muted"], font=FONT,
                           anchor="w")
            lab.grid(row=i, column=0, sticky="w", pady=1)
            self._bind_text(lab, "text", key)
            val = tk.Label(grid, text="—", bg=C["card"], fg=C["text"], font=FONT,
                           anchor="e")
            val.grid(row=i, column=1, sticky="e", pady=1)
            self.est_labels[key] = val

        self.lb_warn = tk.Label(body, text="", bg=C["card"], fg=C["warn"],
                                font=("Segoe UI", 8), anchor="w", justify="left",
                                wraplength=390)
        self.lb_warn.pack(fill="x", pady=(6, 0))

    # ------------------------------------------------------------ 资源条
    def _build_resource_strip(self, root):
        """顶部常驻资源条：CPU / 内存 / 项目盘 / 流程占用，四格横排。"""
        outer = tk.Frame(root, bg=C["border"])
        outer.pack(fill="x", padx=10, pady=(0, 8))
        inner = tk.Frame(outer, bg=C["card"])
        inner.pack(fill="x", padx=1, pady=1)
        self.res = {}
        for i, key in enumerate(("res.cpu", "res.mem", "res.disk", "res.rss")):
            inner.grid_columnconfigure(i, weight=1, uniform="res")
            cell = tk.Frame(inner, bg=C["card"])
            cell.grid(row=0, column=i, sticky="ew",
                      padx=(12 if i == 0 else 10, 10), pady=5)
            lab = tk.Label(cell, text="", bg=C["card"], fg=C["muted"], font=FONT)
            lab.pack(anchor="w")
            self._bind_text(lab, "text", key)
            val = tk.Label(cell, text="n/a", bg=C["card"], fg=C["text"],
                           font=FONT_B, anchor="w")
            val.pack(anchor="w", pady=(0, 2))
            bar = ttk.Progressbar(cell, style="Slim.Horizontal.TProgressbar",
                                  maximum=100, value=0)
            bar.pack(fill="x")
            self.res[key] = (bar, val)

    # ------------------------------------------------------------ 流程卡
    def _build_pipeline(self, parent):
        _o, body, head = self._card(parent, "sec.pipeline", pady=(0, 0),
                                    fill="x", expand=False)
        top = tk.Frame(body, bg=C["card"])
        top.pack(fill="x")
        lab = tk.Label(top, text="", bg=C["card"], fg=C["muted"], font=FONT)
        lab.pack(side="left")
        self._bind_text(lab, "text", "pipe.overall")

        self.lb_pct = tk.Label(top, textvar=self.v_overall, bg=C["card"],
                               fg=C["text"], font=("Segoe UI", 11, "bold"))
        self.lb_pct.pack(side="left", padx=(8, 0))

        # 右侧两组「标签 + 数值」，从右往左：ETA、Elapsed
        for var, key in ((self.v_eta, "pipe.eta"),
                         (self.v_elapsed, "pipe.elapsed")):
            grp = tk.Frame(top, bg=C["card"])
            grp.pack(side="right", padx=(16, 0))
            lab = tk.Label(grp, text="", bg=C["card"], fg=C["muted"], font=FONT)
            lab.pack(side="left", padx=(0, 5))
            self._bind_text(lab, "text", key)
            tk.Label(grp, textvar=var, bg=C["card"], fg=C["text"],
                     font=FONT_B).pack(side="left")

        self.bar = ttk.Progressbar(body, style="Bar.Horizontal.TProgressbar",
                                   maximum=1000, value=0)
        self.bar.pack(fill="x", pady=(6, 4))

        cur_row = tk.Frame(body, bg=C["card"])
        cur_row.pack(fill="x")
        self.lb_current = tk.Label(cur_row, text="", bg=C["card"], fg=C["muted"],
                                   font=FONT, anchor="w")
        self.lb_current.pack(side="left")
        self.lb_meta = tk.Label(cur_row, textvar=self.v_meta, bg=C["card"],
                                fg=C["muted"], font=FONT, anchor="e")
        self.lb_meta.pack(side="right")

        self.board = tk.Canvas(body, bg=C["card"], height=112,
                               highlightthickness=0)
        self.board.pack(fill="x", pady=(6, 4))
        self.board.bind("<Configure>", lambda e: self._render_board())
        self.board.bind("<Button-1>", self._on_board_click)
        self.board.bind("<Motion>", self._board_hover)
        self.board.bind("<Leave>", lambda e: self._hide_tip())
        self._board_ids = {}
        self._hover_sid = None
        self._tip = None
        self._tip_x = self._tip_y = 0

        cols = ("status", "progress", "elapsed", "eta", "check")
        self.tree = ttk.Treeview(body, columns=cols, show="tree headings",
                                 height=5, selectmode="browse")
        self.tree.column("#0", width=190, minwidth=120, stretch=False)
        self.tree.column("status", width=90, anchor="w", stretch=False)
        self.tree.column("progress", width=110, anchor="e", stretch=False)
        self.tree.column("elapsed", width=80, anchor="e", stretch=False)
        self.tree.column("eta", width=80, anchor="e", stretch=False)
        self.tree.column("check", width=110, anchor="w")
        self.tree.pack(fill="x")
        self.tree.tag_configure("running", foreground=C["accent"])
        self.tree.tag_configure("fail", foreground=C["fail"])
        self.tree.tag_configure("check-fail", foreground=C["fail"])
        self.tree.tag_configure("ok", foreground=C["text"])
        self.tree.tag_configure("pending", foreground=C["muted"])
        self.tree.tag_configure("odd", background=C["zebra"])
        self.tree.tag_configure("running", background=C["chip_run"])
        self.tree.tag_configure("fail", background=C["chip_bad"])
        self.tree.tag_configure("check-fail", background=C["chip_bad"])
        self.tree.bind("<<TreeviewSelect>>", self._on_tree_select)

    # ------------------------------------------------------------ 数据库卡
    def _open_db_dialog(self):
        """数据库面板：本地层 + 待下载库 + 大小/状态，独立窗口给足空间。"""
        win = getattr(self, "_db_win", None)
        if win is not None and win.winfo_exists():
            win.lift()
            win.focus_set()
            return
        win = tk.Toplevel(self)
        win.geometry("680x470+%d+%d" % (self.winfo_rootx() + 90,
                                        self.winfo_rooty() + 60))
        win.configure(bg=C["card"])
        win.transient(self)
        self._db_win = win
        self._bind_text(win, "title", "sec.database")

        head = tk.Frame(win, bg=C["card"])
        head.pack(fill="x", padx=12, pady=(10, 6))
        lab = tk.Label(head, text="", bg=C["card"], fg=C["text"], font=FONT_H2)
        lab.pack(side="left")
        self._bind_text(lab, "text", "sec.database")
        b = ttk.Button(head, text="", width=18, command=self._probe_db)
        b.pack(side="right")
        self._bind_text(b, "text", "db.probe")
        self.btn_probe = b
        br = ttk.Button(head, text="", width=9,
                        command=lambda: (self._refresh_db(),
                                         self.db_tree.yview_moveto(0)))
        br.pack(side="right", padx=(0, 6))
        self._bind_text(br, "text", "db.refresh")

        wrap = tk.Frame(win, bg=C["border"])
        wrap.pack(fill="both", expand=True, padx=12)
        self.db_tree = ttk.Treeview(wrap, columns=("size", "id", "state"),
                                    show="tree headings", height=12)
        self.db_tree.column("#0", width=270, minwidth=160)
        self.db_tree.column("size", width=90, anchor="e", stretch=False)
        self.db_tree.column("id", width=110, anchor="w", stretch=False)
        self.db_tree.column("state", width=100, anchor="w", stretch=False)
        self.db_tree.bind("<Double-1>", self._db_open_row)
        self.db_tree.pack(side="left", fill="both", expand=True, padx=1, pady=1)
        sb = ttk.Scrollbar(wrap, orient="vertical", command=self.db_tree.yview)
        sb.pack(side="right", fill="y")
        self.db_tree.configure(yscrollcommand=sb.set)

        foot = tk.Frame(win, bg=C["card"])
        foot.pack(fill="x", padx=12, pady=(6, 2))
        self.lb_dbtotal = tk.Label(foot, textvar=self.v_db_total, bg=C["card"],
                                   fg=C["text"], font=FONT_B, anchor="w")
        self.lb_dbtotal.pack(side="left")
        cb = ttk.Button(foot, text="", width=8, command=win.destroy)
        cb.pack(side="right")
        self._bind_text(cb, "text", "db.close")
        bc = ttk.Button(foot, text="", width=11, command=self._copy_db_path)
        bc.pack(side="right", padx=(0, 6))
        self._bind_text(bc, "text", "db.copy_path")

        hint = tk.Label(win, text="", bg=C["card"], fg=C["hint"],
                        font=("Segoe UI", 8), anchor="w", justify="left",
                        wraplength=640)
        hint.pack(fill="x", padx=12, pady=(0, 10))
        self._bind_text(hint, "text", "db.hint")

        self._render_db()
        self._refresh_db()
        self.db_tree.yview_moveto(0)

    # ------------------------------------------------------------ 日志卡
    def _build_log(self, parent):
        head = tk.Frame(parent, bg=C["card"])
        head.pack(fill="x", padx=10, pady=(8, 4))
        body = tk.Frame(parent, bg=C["card"])
        body.pack(fill="both", expand=True, padx=10, pady=(0, 8))

        # 查找条（Ctrl+F 显示）
        self.find_bar = tk.Frame(body, bg=C["card"])
        inb = tk.Frame(self.find_bar, bg=C["card"])
        inb.pack(fill="x")
        lab = tk.Label(inb, text="", bg=C["card"], fg=C["muted"], font=FONT)
        lab.pack(side="left", padx=(0, 4))
        self._bind_text(lab, "text", "log.find")
        self.v_find = tk.StringVar()
        self.ent_find = ttk.Entry(inb, textvariable=self.v_find, width=22)
        self.ent_find.pack(side="left")
        self.ent_find.bind("<KeyRelease>", self._find_refresh)
        self.ent_find.bind("<Return>", lambda e: self._find_step(1))
        self.ent_find.bind("<Shift-Return>", lambda e: self._find_step(-1))
        self.ent_find.bind("<Escape>", lambda e: self._toggle_find(False))
        self.lb_find = tk.Label(inb, text="0/0", bg=C["card"], fg=C["muted"],
                                font=FONT, width=7)
        self.lb_find.pack(side="left", padx=6)
        for key, step in (("log.find_prev", -1), ("log.find_next", 1)):
            fb = ttk.Button(inb, text="", width=5,
                            command=lambda s=step: self._find_step(s))
            fb.pack(side="left", padx=(0, 4))
            self._bind_text(fb, "text", key)
        bx = ttk.Button(inb, text="✕", width=3,
                        command=lambda: self._toggle_find(False))
        bx.pack(side="left")
        self._find_hits = []
        self._find_i = -1

        ttl = tk.Label(head, text="", bg=C["card"], fg=C["text"], font=FONT_H2)
        ttl.pack(side="left", padx=(0, 14))
        self._bind_text(ttl, "text", "sec.log")
        lab = tk.Label(head, text="", bg=C["card"], fg=C["muted"], font=FONT)
        lab.pack(side="left", padx=(0, 4))
        self._bind_text(lab, "text", "log.level")
        self.cb_level = ttk.Combobox(head, width=12, state="readonly")
        self.cb_level.pack(side="left", padx=(0, 10))
        self.cb_level.bind("<<ComboboxSelected>>", self._on_level)
        self._level_label_keys = ("log.all", "log.warn_plus", "log.errors")
        self._render_level_values()

        lab = tk.Label(head, text="", bg=C["card"], fg=C["muted"], font=FONT)
        lab.pack(side="left", padx=(0, 4))
        self._bind_text(lab, "text", "log.filter")
        ent = ttk.Entry(head, textvariable=self.v_filter, width=18)
        ent.pack(side="left")
        ent.bind("<KeyRelease>", lambda e: self._render_log(force=True))
        ent.bind("<FocusOut>", lambda e: self._save_cfg())
        ent.bind("<Return>", lambda e: self._save_cfg())
        self.ent_filter = ent

        cb = tk.Checkbutton(head, text="", variable=self.v_autoscroll,
                            bg=C["card"], fg=C["muted"], font=FONT, bd=0,
                            selectcolor=C["field"], highlightthickness=0,
                            activebackground=C["card"])
        cb.pack(side="right")
        self._bind_text(cb, "text", "log.autoscroll")
        b = ttk.Button(head, text="", width=8,
                       command=lambda: (self.txt.delete("1.0", "end"),
                                        setattr(self, "_log_lines", 0)))
        b.pack(side="right", padx=(0, 6))
        self._bind_text(b, "text", "log.clear")

        wrap = tk.Frame(body, bg=C["border"])
        wrap.pack(fill="both", expand=True)
        self._log_wrap = wrap
        self.txt = tk.Text(wrap, font=FONT_MONO, bd=0, wrap="word",
                           bg=C["log_bg"], fg=C["log_fg"], insertwidth=0,
                           highlightthickness=0, state="disabled")
        self.txt.pack(side="left", fill="both", expand=True, padx=1, pady=1)
        sb = ttk.Scrollbar(wrap, orient="vertical", command=self.txt.yview)
        sb.pack(side="right", fill="y")
        self.txt.configure(yscrollcommand=sb.set)
        # 用户往上翻时自动关掉"自动滚动"，翻回底部再打开
        self.txt.bind("<MouseWheel>", self._log_scrolled, add="+")
        self.txt.bind("<ButtonRelease-1>", self._log_scrolled, add="+")
        self.txt.bind("<KeyRelease>", self._log_scrolled, add="+")
        sb.bind("<ButtonRelease-1>", self._log_scrolled, add="+")
        self.txt.bind("<Button-3>", self._log_menu)
        self._menu = tk.Menu(self, tearoff=0, bg=C["card"], fg=C["text"],
                             activebackground=C["sel"],
                             activeforeground=C["text"], bd=1,
                             relief="solid")
        self._menu.add_command(label="", command=self._copy_selection)
        self._menu.add_command(label="", command=self._copy_all)
        self._menu.add_separator()
        self._menu.add_command(label="", command=self._filter_from_selection)
        self._menu.add_command(label="", command=self._open_log_file)
        self._render_log_menu_labels()
        self.txt.tag_configure("warn", foreground=C["warn"])
        self.txt.tag_configure("err", foreground=C["fail"])
        self.txt.tag_configure("ok", foreground=C["ok"])
        self.txt.tag_configure("dim", foreground=C["hint"])
        self.txt.tag_configure("hit", background=C["hit"])
        self.txt.tag_configure("find", background=C["find_bg"])
        self.txt.tag_configure("find_cur", background=C["accent"],
                               foreground=C["field"])
        self._log_data = []

    def _render_level_values(self):
        if not hasattr(self, "cb_level"):
            return
        vals = [self.t(k) for k in self._level_label_keys]
        self.cb_level.configure(values=vals)
        self.cb_level.current(min(self._level_idx, len(vals) - 1))

    def _on_font_scale(self, _ev=None):
        try:
            pct = int(str(self.cb_font.get()).strip().rstrip("%"))
        except ValueError:
            return
        self._apply_font_scale(pct / 100.0)

    def _on_level(self, _ev=None):
        idx = self.cb_level.current()
        self._level_idx = max(0, min(2, idx))
        self._render_log(force=True)
        self._save_cfg()

    # ------------------------------------------------------------ 状态栏
    def _build_statusbar(self, root):
        bar = tk.Frame(root, bg=C["card"], height=40)
        bar.pack(fill="x", side="bottom", padx=10, pady=(0, 6))
        bar.pack_propagate(False)

        self.btn_start = ttk.Button(bar, text="", width=14,
                                    command=lambda: self._start(None))
        self.btn_start.pack(side="left", padx=(12, 6), pady=8)
        self._bind_text(self.btn_start, "text", "btn.start")
        self.btn_resume = ttk.Button(bar, text="", width=12,
                                     command=lambda: self._start("resume"))
        self.btn_resume.pack(side="left", padx=(0, 6), pady=8)
        self._bind_text(self.btn_resume, "text", "btn.resume")
        self.btn_stop = ttk.Button(bar, text="", width=10, command=self._stop,
                                   state="disabled")
        self.btn_stop.pack(side="left", padx=(0, 6), pady=8)
        self._bind_text(self.btn_stop, "text", "btn.stop")
        b = ttk.Button(bar, text="", width=16, command=self._open_results)
        b.pack(side="left", padx=(0, 6), pady=8)
        self._bind_text(b, "text", "btn.results")
        b2 = ttk.Button(bar, text="", width=16, command=self._export_summary)
        b2.pack(side="left", pady=8)
        self._bind_text(b2, "text", "btn.export")

        self.lb_status = tk.Label(bar, text="", bg=C["card"], fg=C["muted"],
                                  font=FONT, anchor="e")
        self.lb_status.pack(side="right", padx=14)

    # ------------------------------------------------------------ 数据回流
    def _restore(self, prefill_workdir):
        for p in self.cfg.get("raws") or []:
            if os.path.exists(p):
                self.v_raws.append(p)
        if prefill_workdir:
            self.v_workdir.set(os.path.abspath(prefill_workdir))
        self._refresh_rawsum()
        self._refresh_db()
        wd = self.v_workdir.get().strip()
        if wd and os.path.isdir(wd):
            self._attach(wd, quiet=True)

    def _refresh_rawsum(self):
        n, total = model.raw_stats(self.v_raws)
        if n and total:
            self.v_rawsum.set(self.t("inp.summary", n=n, size=fmt_bytes(total)))
        elif n:
            self.v_rawsum.set(self.t("inp.summary_n", n=n))
        else:
            self.v_rawsum.set(self.t("inp.no_raw"))

    def _render_rawlist(self):
        self.lb_raw.delete(0, "end")
        for p in self.v_raws:
            self.lb_raw.insert("end", os.path.basename(p) or p)
        self._refresh_rawsum()

    def _add_files(self):
        paths = filedialog.askopenfilenames(
            title=self.t("inp.add_files"),
            filetypes=[("Thermo RAW", "*.raw"), ("All files", "*.*")])
        self._add(paths)

    def _add_folder(self):
        d = filedialog.askdirectory(title=self.t("inp.add_folder"))
        if d:
            self._add([d])

    def _add(self, paths):
        for p in paths or []:
            p = os.path.abspath(p)
            if p not in self.v_raws:
                self.v_raws.append(p)
        self._render_rawlist()
        self._refresh_estimate()
        self._save_cfg()

    def _remove_selected(self):
        for i in sorted(self.lb_raw.curselection(), reverse=True):
            del self.v_raws[i]
        self._render_rawlist()
        self._refresh_estimate()
        self._save_cfg()

    def _clear_raws(self):
        self.v_raws = []
        self._render_rawlist()
        self._refresh_estimate()
        self._save_cfg()

    def _pick_workdir(self):
        d = filedialog.askdirectory(title=self.t("inp.workdir"))
        if d:
            self.v_workdir.set(os.path.abspath(d))
            self._refresh_db()
            self._attach(d, quiet=True)
            self._refresh_estimate()
            self._save_cfg()

    def _load_existing(self):
        d = filedialog.askdirectory(title=self.t("inp.attach"))
        if not d:
            return
        self.v_workdir.set(os.path.abspath(d))
        ok = self._attach(d, quiet=False)
        if ok:
            self._refresh_db()
            self._refresh_estimate()

    def _attach(self, workdir, quiet=True):
        """载入（或接管）项目里的 progress.jsonl。"""
        path = os.path.join(workdir, "00_logs", "progress.jsonl")
        if not os.path.isfile(path):
            if not quiet:
                messagebox.showinfo("ORFscape", self.t("inp.notfound"))
            return False
        self._events_path = path
        self._ev_off = 0
        self.state = runlog.RunState()
        recs, self._ev_off = runlog.read_events(path, 0)
        for r in recs:
            self.state.feed(r)
        self._log_data = []
        self.txt.configure(state="normal")
        self.txt.delete("1.0", "end")
        self.txt.configure(state="disabled")
        self._log_lines = 0
        for r in recs:
            if r.get("ev") == "log":
                self._append_log(r.get("line") or "", r.get("stage"))
        self._apply(recs, quiet=quiet)
        if recs and not quiet:
            self._set_status(self.t("inp.loaded", n=len(recs),
                                    path=os.path.basename(os.path.dirname(path))))
        return True

    def _save_cfg(self):
        self.cfg.update({
            "raws": list(self.v_raws),
            "species": self.v_species.get().strip(),
            "desc": self.v_desc.get().strip(),
            "workdir": self.v_workdir.get().strip(),
            "threads": int(self.v_threads.get() or 16),
            "workers": int(self.v_workers.get() or 4),
            "dbs": self.v_dbs.get().strip(),
            "skip_p6b": bool(self.v_skip_p6b.get()),
            "skip_qc": bool(self.v_skip_qc.get()),
            "force": bool(self.v_force.get()),
            "lang": self.lang,
            "advanced": bool(self._adv_open.get()),
            "theme": theme.get_theme(),
            "confirm": bool(self.v_confirm.get()),
            "filter": self.v_filter.get(),
            "log_level": int(getattr(self, "_level_idx", 0)),
            "font_scale": float(getattr(self, "_font_scale", 1.0)),
        })
        model.save_settings(self.cfg)

    def _on_lang(self, _ev=None):
        idx = self.cb_lang.current()
        self.lang = i18n.LANGS[max(0, idx)][0]
        self.retranslate()
        self._save_cfg()

    # ------------------------------------------------------------ 估算
    def _form(self):
        return {
            "raws": list(self.v_raws),
            "species": self.v_species.get().strip(),
            "desc": self.v_desc.get().strip(),
            "threads": int(self.v_threads.get() or 16),
            "workers": int(self.v_workers.get() or 4),
            "dbs": self.v_dbs.get().strip(),
            "skip_p6b": bool(self.v_skip_p6b.get()),
            "skip_qc": bool(self.v_skip_qc.get()),
            "force": bool(self.v_force.get()),
        }

    def _refresh_estimate(self):
        form = self._form()
        wd = self.v_workdir.get().strip()

        def work():
            try:
                res = estimator.estimate(wd, form["raws"], threads=form["threads"],
                                         workers=form["workers"], dbs=form["dbs"])
            except Exception as exc:              # 估算失败不该卡界面
                res = {"error": str(exc)}
            self._ui_q.put((self._on_estimate, res))

        threading.Thread(target=work, daemon=True).start()

    def _on_estimate(self, res):
        self._est = res
        if res.get("error"):
            self._set_est_text({})
            return
        i = res["inputs"]
        g = res["genome"]
        t = res["tasks"]
        d = res["disk"]
        self._set_est_text({
            "est.raw": "%d · %s" % (i["n_raw"], fmt_bytes(i["bytes"])),
            "est.ms2": self.t("est.value_src", v="{:,}".format(i["ms2"]),
                              src=self.t("est.src." + i["source"])),
            "est.genome": self.t("est.value_src", v="%.2f Mb" % g["mb"],
                                 src=self.t("est.src." + g["source"])),
            "est.tasks": "%d / %d" % (t["e1"], t["e2"]),
            "est.time": self.t("est.value_src", v=fmt_sec(res["total_sec"]),
                               src=self.t("est.src.calibrated" if res["calibrated"]
                                          else "est.src.estimated")),
            "est.disk": self.t("est.disk_val", need=fmt_bytes(d["total"]),
                               free=fmt_bytes((d.get("free_gb") or 0) * estimator.GB)),
            "est.ram": "%.1f GB" % res["ram"]["need_gb"],
        })
        warns = [self.t(w["key"], **w["args"]) for w in res.get("warnings") or []]
        self.lb_warn.configure(text="\n".join(warns), fg=C["warn"])
        self._sampler.set_disk_path(self.v_workdir.get().strip() or os.getcwd())

    def _set_est_text(self, mapping):
        for key, lab in self.est_labels.items():
            lab.configure(text=mapping.get(key, "—"))

    # ------------------------------------------------------------ 数据库
    def _refresh_db(self):
        wd = self.v_workdir.get().strip()
        sp = self.v_species.get().strip()
        code = sp if sp and " " not in sp else None
        try:
            rows = estimator.database_rows(wd, code, sp or None) if wd else []
        except Exception:
            rows = []
        self._db_rows = rows
        self._render_db()

    def _render_db(self):
        if not hasattr(self, "db_tree"):
            return
        self.db_tree.heading("#0", text=self.t("db.col.name"))
        self.db_tree.heading("size", text=self.t("db.col.size"))
        self.db_tree.heading("id", text=self.t("db.col.id"))
        self.db_tree.heading("state", text=self.t("db.col.state"))
        self.db_tree.delete(*self.db_tree.get_children())
        total = 0
        for r in self._db_rows:
            size = r.get("bytes")
            if size:
                total += size
            ident = r.get("accession") or r.get("proteome") or ""
            self.db_tree.insert(
                "", "end", text=r.get("label"), iid=r.get("key"),
                values=(fmt_bytes(size) if size else "—", ident,
                        self.t("db.state." + (r.get("state") or "unknown"))))
        self.v_db_total.set("%s  %s" % (self.t("db.total"), fmt_bytes(total))
                            if self._db_rows else "")
        try:
            self.db_tree.update_idletasks()
            self.db_tree.yview_moveto(0)
        except tk.TclError:
            pass

    def _probe_db(self):
        wd = self.v_workdir.get().strip()
        if not wd:
            messagebox.showwarning("ORFscape", self.t("msg.need_dir"))
            return
        sp = self.v_species.get().strip() or None
        if hasattr(self, "btn_probe"):
            self.btn_probe.configure(state="disabled", text=self.t("db.probing"))

        def work():
            try:
                rows = estimator.probe_remote_sizes(self._db_rows, sp, log=lambda *_: None)
            except Exception:
                rows = self._db_rows
            self._ui_q.put((self._on_probe, rows))

        threading.Thread(target=work, daemon=True).start()

    def _on_probe(self, rows):
        self._db_rows = rows or self._db_rows
        if hasattr(self, "btn_probe"):
            self.btn_probe.configure(state="normal")
            self.btn_probe.configure(text=self.t("db.probe"))
        self._render_db()

    def _db_row(self, key):
        for r in self._db_rows:
            if r.get("key") == key:
                return r
        return None

    def _db_selected(self):
        sel = self.db_tree.selection() if hasattr(self, "db_tree") else ()
        return self._db_row(sel[0]) if sel else None

    def _copy_db_path(self):
        row = self._db_selected()
        if not row:
            return
        text = row.get("url") or row.get("path") or row.get("label") or ""
        if not text:
            return
        self.clipboard_clear()
        self.clipboard_append(str(text))
        self._set_status(self.t("msg.copied", what=str(text)), color=C["ok"])

    def _db_open_row(self, _ev=None):
        """双击：本地文件打开所在目录；远程项退化为复制链接。"""
        row = self._db_selected()
        if not row:
            return
        path = row.get("path")
        if path and os.path.isfile(path):
            target = os.path.dirname(path)
            try:
                if os.name == "nt":
                    os.startfile(target)             # noqa: S606
                elif sys.platform == "darwin":
                    subprocess.Popen(["open", target])
                else:
                    subprocess.Popen(["xdg-open", target])
            except OSError:
                pass
        else:
            self._copy_db_path()

    # ------------------------------------------------------------ 画板
    def _render_board(self):
        if not hasattr(self, "board"):
            return
        cv = self.board
        cv.delete("all")
        self._board_ids = {}
        self._hover_sid = None
        self._tip = None
        self._tip_x = self._tip_y = 0
        w = max(320, cv.winfo_width() or 900)
        h = int(cv.cget("height"))
        order = self.state.order or runlog.DEFAULT_ORDER
        n = len(order)
        per_row = 7 if n > 7 else n
        pad_x, pad_y, gap_x, gap_y = 2, 4, 9, 14
        rows = (n + per_row - 1) // per_row
        ch = int((h - 2 * pad_y - (rows - 1) * gap_y) / max(1, rows))
        cw = min((w - 2 * pad_x - (per_row - 1) * gap_x) / float(per_row), 152.0)

        for idx, sid in enumerate(order):
            r, c = divmod(idx, per_row)
            in_row = min(per_row, n - r * per_row)
            row_w = in_row * cw + (in_row - 1) * gap_x
            x_off = (w - row_w) / 2.0
            x0 = x_off + c * (cw + gap_x)
            y0 = pad_y + r * (ch + gap_y)
            x1, y1 = x0 + cw, y0 + ch
            st = self.state.stages.get(sid)
            status = st.status if st else "pending"
            border, fill, dot = st_style(status)

            # 行内连接：短线 + 小箭头，标出执行方向
            if c < in_row - 1:
                cy = y0 + ch / 2.0
                cv.create_line(x1 + 1, cy, x1 + gap_x - 3, cy,
                               fill=C["border"], width=1)
                cv.create_polygon(x1 + gap_x - 4, cy - 2.5,
                                  x1 + gap_x - 4, cy + 2.5,
                                  x1 + gap_x - 1, cy,
                                  fill=C["border"], outline="")

            _rr(cv, x0, y0, x1, y1, 7, fill=fill, outline=border,
                width=2 if status == "running" else 1)
            cv.create_text(x0 + 9, y0 + 13, text=sid, anchor="w",
                           font=self._fnt(9, "bold"),
                           fill=C["accent"] if status == "running" else C["text"])
            title = self.t("stage_short." + sid)
            cv.create_text(x0 + 9, y0 + ch / 2 + 3, text=title, anchor="w",
                           font=self._fnt(8), fill=C["muted"])

            # 右上角状态：勾/叉/点
            gx, gy = x1 - 11, y0 + 12
            if status == "ok":
                cv.create_text(gx, gy, text="✓", font=self._fnt(9, "bold"),
                               fill=C["ok"])
            elif status in ("fail", "check-fail"):
                cv.create_text(gx, gy, text="✕", font=self._fnt(9, "bold"),
                               fill=C["fail"])
            elif status == "running":
                cv.create_oval(gx - 4, gy - 4, gx + 4, gy + 4, fill=dot,
                               outline=dot)
            elif status in ("skip", "aborted"):
                cv.create_line(gx - 4, gy, gx + 4, gy, fill=C["idle"], width=2)
            else:
                cv.create_oval(gx - 3, gy - 3, gx + 3, gy + 3, fill="",
                               outline=C["border"], width=1)

            frac = st.fraction() if st else 0.0
            bx0, bx1, by = x0 + 8, x1 - 8, y1 - 5
            cv.create_rectangle(bx0, by, bx1, by + 2, fill=C["track"], outline="")
            if frac > 0:
                cv.create_rectangle(bx0, by, bx0 + (bx1 - bx0) * frac, by + 2,
                                    fill=C["accent"] if status == "running" else dot,
                                    outline="")
            if status == "running" and st and st.total:
                cv.create_text(x1 - 8, y1 - 10,
                               text="%d%%" % int(100 * frac), anchor="se",
                               font=self._fnt(8, "bold"), fill=C["accent"])
            self._board_ids[sid] = (x0, y0, x1, y1)

    def _show_tip(self, sid):
        st = self.state.stages.get(sid)
        if st is None:
            return
        est = ((self._est or {}).get("est_sec", {}) or {}).get(sid) or {}
        rows = model.stage_tooltip(st, est.get("sec"))
        if not rows:
            return
        if self._tip is None or not self._tip.winfo_exists():
            self._tip = tk.Toplevel(self)
            self._tip.overrideredirect(True)
            self._tip.configure(bg=C["border"])
            self.lb_tip = tk.Label(self._tip, text="", bg=C["card"],
                                   fg=C["text"], font=FONT, justify="left",
                                   anchor="w", padx=8, pady=6)
            self.lb_tip.pack(padx=1, pady=1)
        title = "%s · %s" % (sid, self.t("stage." + sid))
        lines = [title]
        for key, args in rows:
            lines.append(self.t(key, **args))
        self.lb_tip.configure(text=chr(10).join(lines), bg=C["card"],
                              fg=C["text"])
        x = self.board.winfo_rootx() + self._tip_x + 14
        y = self.board.winfo_rooty() + self._tip_y + 18
        self._tip.geometry("+%d+%d" % (x, y))
        self._tip.deiconify()

    def _hide_tip(self):
        self._hover_sid = None
        if self._tip is not None:
            try:
                self._tip.destroy()
            except tk.TclError:
                pass
            self._tip = None

    def _board_hover(self, ev):
        sid = None
        for key, (x0, y0, x1, y1) in getattr(self, "_board_ids", {}).items():
            if x0 <= ev.x <= x1 and y0 <= ev.y <= y1:
                sid = key
                break
        self._tip_x, self._tip_y = ev.x, ev.y
        if sid == self._hover_sid:
            return
        self._hover_sid = sid
        if not sid:
            self._hide_tip()
            self._hover_sid = None
            return
        self._show_tip(sid)

    def _on_board_click(self, ev):
        for sid, (x0, y0, x1, y1) in getattr(self, "_board_ids", {}).items():
            if x0 <= ev.x <= x1 and y0 <= ev.y <= y1:
                self.v_filter.set(sid)
                self._render_log(force=True)
                return

    # ------------------------------------------------------------ 表格
    def _render_table_headings(self):
        if not hasattr(self, "tree"):
            return
        self.tree.heading("#0", text=self.t("pipe.col.stage"))
        for col in ("status", "progress", "elapsed", "eta", "check"):
            self.tree.heading(col, text=self.t("pipe.col." + col))

    def _render_table(self):
        """增量刷新：只改值变了的行，不整表重建（每 400 ms 调一次）。"""
        if not hasattr(self, "tree"):
            return
        self._render_table_headings()
        rows = self.state.rows()
        ids = [st.sid for st in rows]
        existing = list(self.tree.get_children())
        cache = getattr(self, "_table_cache", None)
        if cache is None or ids != existing:
            self.tree.delete(*existing)
            self._table_cache = {}
            cache = self._table_cache
            for idx, st in enumerate(rows):
                tags = (st.status, "odd") if idx % 2 else (st.status,)
                self.tree.insert("", "end", iid=st.sid, text="", values=(),
                                 tags=tags)
        for st in rows:
            prog = ""
            if st.total:
                prog = "%d/%d" % (st.done, st.total)
            elif st.status == "running":
                prog = "…"
            check = ""
            if st.check:
                check = "%dF/%dW" % (st.check.get("fails") or 0,
                                     st.check.get("warns") or 0)
            rem = st.remaining_sec()
            text = "%s  %s" % (st.sid, self.t("stage." + st.sid))
            values = (self.t("st." + st.status), prog,
                      fmt_sec(st.elapsed) if st.elapsed else "",
                      fmt_sec(rem) if (st.status == "running" and rem) else "",
                      check)
            sig = (text, values)
            if cache.get(st.sid) == sig:
                continue
            cache[st.sid] = sig
            if self.tree.exists(st.sid):
                self.tree.item(st.sid, text=text, values=values)
        cur = self.state.failed() or self.state.current()
        if cur is not None and self.tree.exists(cur.sid)                 and getattr(self, "_table_seen", None) != cur.sid:
            self._table_seen = cur.sid
            self.tree.see(cur.sid)

    def _on_tree_select(self, _ev=None):
        sel = self.tree.selection()
        if sel:
            self.v_filter.set(sel[0])
            self._render_log(force=True)

    # ------------------------------------------------------------ 刷新
    def _apply(self, recs, quiet=False):
        st = self.state
        frac = st.overall()
        self._set_bar(int(frac * 1000))
        self.v_overall.set("%.0f%%" % (frac * 100))
        self._set_title_progress(frac if st.status == "running" else None)
        if st.started:
            el = (st.ended or datetime.datetime.now().timestamp()) - st.started
            self.v_elapsed.set(fmt_sec(el))
        eta = st.eta_sec()
        self.v_eta.set(self.t("pipe.eta_unknown") if eta is None else fmt_sec(eta))
        cur = st.current()
        if st.status in runlog.TERMINAL:
            self.v_current.set(self._status_text(st.status))
        elif cur is not None:
            extra = ""
            if cur.total:
                extra = " · " + self.t("pipe.step", done=cur.done,
                                       total=cur.total, unit=cur.unit)
            elif cur.psm:
                extra = " · " + self.t("pipe.psm", n=cur.psm)
            self.v_current.set("%s · %s%s" % (self.t("pipe.current"),
                                              self.t("stage." + cur.sid), extra))
        else:
            key = "pipe.idle" if (st.n_events or self.v_workdir.get().strip())                 else "pipe.empty"
            self.v_current.set(self.t(key))
        done = sum(1 for x in st.rows() if x.status in ("ok", "skip"))
        psm = st.psm_total()
        bits = [self.t("pipe.stages", done=done, total=len(st.order))]
        if psm:
            bits.append(self.t("pipe.psm", n="{:,}".format(psm)))
        self.v_meta.set("  ·  ".join(bits))
        self._render_board()
        self._render_table()
        if quiet:
            return

    def _status_text(self, status):
        st = self.state
        el = fmt_sec((st.ended or 0) - (st.started or 0)) if st.started else "—"
        cur = st.failed() or st.current()
        name = self.t("stage." + cur.sid) if cur else "—"
        if status == "ok":
            return self.t("msg.finished", t=el)
        if status == "fail":
            return self.t("msg.failed", stage=name)
        if status == "check-fail":
            return self.t("msg.checkfail", stage=name)
        if status == "stopped":
            return self.t("msg.stopped", t=el)
        return ""

    # ------------------------------------------------------------ 日志
    def _append_log(self, line, stage=None):
        line = (line or "").rstrip("\n")
        if not line:
            return
        self._log_data.append((line, stage))
        if len(self._log_data) > 6000:
            self._log_data = self._log_data[-5000:]
        self._render_log()

    @staticmethod
    def _line_level(line):
        return model.line_level(line)

    def _render_log(self, force=False):
        if not hasattr(self, "txt"):
            return
        nl = chr(10)
        pat = (self.v_filter.get() or "").strip().lower()
        by_stage = pat in {x.lower() for x in self.state.order}
        lines = model.filter_log(self._log_data, pat,
                                 getattr(self, "_level_idx", 0),
                                 self.state.order)
        self.txt.configure(state="normal")
        self.txt.delete("1.0", "end")
        for line, lv in lines[-3000:]:
            tag = {2: "err", 1: "warn", 0: ""}[lv]
            if not tag and (" ok" in line.lower() or "done" in line.lower()
                            or "完成" in line):
                tag = "ok"
            start = self.txt.index("end-1c")
            self.txt.insert("end", line + nl, tag)
            if pat and not by_stage:
                pos = start
                for _ in range(40):          # 单行最多 40 处高亮
                    pos = self.txt.search(pat, pos, stopindex=start + " lineend",
                                          nocase=True)
                    if not pos:
                        break
                    endp = "%s+%dc" % (pos, len(pat))
                    self.txt.tag_add("hit", pos, endp)
                    pos = endp
        self.txt.insert("end", nl)       # 末尾留一行，避免最后一行贴边
        self.txt.configure(state="disabled")
        self._log_lines = len(lines)
        if self.v_autoscroll.get():
            self.txt.see("end")


    # ------------------------------------------------------------ 日志：查找与右键
    def _toggle_find(self, show=None):
        """Ctrl+F 显示/隐藏查找条。"""
        if not hasattr(self, "find_bar"):
            return "break"
        showing = bool(self.find_bar.winfo_manager())
        if show is None:
            show = not showing
        if show:
            self.find_bar.pack(fill="x", pady=(0, 4), before=self._log_wrap)
            self.ent_find.focus_set()
            self.ent_find.select_range(0, "end")
            self._find_refresh()
        else:
            self.find_bar.pack_forget()
            for tag in ("find", "find_cur"):
                self.txt.tag_remove(tag, "1.0", "end")
            self._find_hits = []
            self._find_i = -1
        return "break"

    def _find_refresh(self, _ev=None):
        needle = self.v_find.get()
        for tag in ("find", "find_cur"):
            self.txt.tag_remove(tag, "1.0", "end")
        hits = []
        if needle:
            pos = "1.0"
            for _ in range(5000):                 # 上限，避免超长日志卡死
                pos = self.txt.search(needle, pos, stopindex="end", nocase=True)
                if not pos:
                    break
                end = "%s+%dc" % (pos, len(needle))
                hits.append(pos)
                self.txt.tag_add("find", pos, end)
                pos = end
        self._find_hits = hits
        self._find_i = -1
        self._find_step(1)

    def _find_step(self, step):
        hits = getattr(self, "_find_hits", [])
        if not hits:
            if hasattr(self, "lb_find"):
                self.lb_find.configure(text="0/0")
            return
        self._find_i = (self._find_i + step) % len(hits)
        pos = hits[self._find_i]
        end = "%s+%dc" % (pos, len(self.v_find.get()))
        self.txt.tag_remove("find_cur", "1.0", "end")
        self.txt.tag_add("find_cur", pos, end)
        self.txt.tag_raise("find_cur")
        self.txt.see(pos)
        self.lb_find.configure(text="%d/%d" % (self._find_i + 1, len(hits)))

    def _render_log_menu_labels(self):
        if not hasattr(self, "_menu"):
            return
        for idx, key in ((0, "menu.copy"), (1, "menu.copy_all"),
                         (3, "menu.filter"), (4, "menu.open_log")):
            try:
                self._menu.entryconfigure(idx, label=self.t(key))
            except tk.TclError:
                pass

    def _log_menu(self, ev):
        try:
            self._menu.tk_popup(ev.x_root, ev.y_root)
        finally:
            try:
                self._menu.grab_release()
            except tk.TclError:
                pass

    def _selection(self):
        try:
            return self.txt.get("sel.first", "sel.last")
        except tk.TclError:
            return ""

    def _copy_selection(self):
        text = self._selection()
        if text:
            self.clipboard_clear()
            self.clipboard_append(text)

    def _copy_all(self):
        self.clipboard_clear()
        self.clipboard_append(chr(10).join(line for line, _st in self._log_data))

    def _filter_from_selection(self):
        text = self._selection().strip().split(chr(10))[0].strip()
        if text:
            self.v_filter.set(text[:120])
            self._render_log(force=True)

    def _open_log_file(self):
        wd = self.v_workdir.get().strip()
        path = model.latest_log(wd) if wd else None
        if not path:
            self._set_status(self.t("msg.no_log"), color=C["warn"])
            return
        try:
            if os.name == "nt":
                os.startfile(str(path))               # noqa: S606
            elif sys.platform == "darwin":
                subprocess.Popen(["open", str(path)])
            else:
                subprocess.Popen(["xdg-open", str(path)])
        except OSError as exc:
            messagebox.showerror("ORFscape", str(exc))

    def _log_scrolled(self, _ev=None):
        """往上翻日志时暂停跟随；翻回底部自动恢复。"""
        if not hasattr(self, "txt"):
            return
        try:
            _top, bottom = self.txt.yview()
        except tk.TclError:
            return
        at_bottom = bottom >= 0.999
        if at_bottom != bool(self.v_autoscroll.get()):
            self.v_autoscroll.set(at_bottom)

    def _set_bar(self, target):
        """总进度条平滑过渡（8 帧左右），避免阶段切换时跳变。"""
        self._bar_target = target
        if not hasattr(self, "bar"):
            return
        try:
            cur = float(self.bar.cget("value"))
        except (tk.TclError, ValueError):
            return
        if abs(cur - target) < 1:
            self.bar.configure(value=target)
            return
        if not getattr(self, "_bar_anim", False):
            self._bar_anim = True
            self.after(30, self._animate_bar)

    def _animate_bar(self):
        if not hasattr(self, "bar"):
            self._bar_anim = False
            return
        try:
            cur = float(self.bar.cget("value"))
        except (tk.TclError, ValueError):
            self._bar_anim = False
            return
        target = getattr(self, "_bar_target", cur)
        if abs(target - cur) <= 1:
            self.bar.configure(value=target)
            self._bar_anim = False
            return
        self.bar.configure(value=cur + (target - cur) / 3.0)
        self.after(30, self._animate_bar)

    def _set_title_progress(self, frac):
        """运行中把百分比写进标题；任务栏进度不可用时它就是降级显示。"""
        if frac is None:
            if self._title_pct is not None:
                self._title_pct = None
                self.title("ORFscape")
            return
        try:
            pct = int(round(max(0.0, min(1.0, float(frac))) * 100))
        except (TypeError, ValueError):
            return
        self._taskbar.set_progress(pct / 100.0)
        if pct == self._title_pct:
            return
        self._title_pct = pct
        self.title("ORFscape — %d%%" % pct)

    def _flash_taskbar(self):
        """Windows：跑完闪任务栏（窗口不在前台时才有意义）。"""
        if os.name != "nt":
            return
        try:
            if self.focus_displayof() is not None:
                return
        except tk.TclError:
            pass
        try:
            import ctypes
            from ctypes import wintypes

            class _FLASHWINFO(ctypes.Structure):
                _fields_ = [("cbSize", wintypes.UINT), ("hwnd", wintypes.HWND),
                            ("dwFlags", wintypes.DWORD), ("uCount", wintypes.UINT),
                            ("dwTimeout", wintypes.DWORD)]

            GA_ROOT = 2
            FLASHW_TRAY = 2
            FLASHW_TIMERNOFG = 12
            hwnd = ctypes.windll.user32.GetAncestor(self.winfo_id(), GA_ROOT)
            info = _FLASHWINFO(ctypes.sizeof(_FLASHWINFO), hwnd,
                               FLASHW_TRAY | FLASHW_TIMERNOFG, 3, 0)
            ctypes.windll.user32.FlashWindowEx(ctypes.byref(info))
        except Exception:
            pass

    # ------------------------------------------------------------ 资源
    def _sample_resources(self):
        wd = self.v_workdir.get().strip()
        if wd and os.path.isdir(wd):
            self._sampler.set_disk_path(wd)
        pid = self._proc.pid if self._proc else None
        self._sample = self._sampler.sample(pid)
        s = self._sample
        cpu = s.get("cpu_pct")
        cpu_txt = (self.t("res.na") if cpu is None else
                   "%d%% · %s" % (int(cpu), self.t("res.cores",
                                                   n=s.get("cpu_n") or 0)))
        self._set_res("res.cpu", cpu, cpu_txt)
        mp = s.get("mem_pct")
        self._set_res("res.mem", mp,
                      self.t("res.mem_val", used=s.get("mem_used_gb"),
                             total=s.get("mem_total_gb"))
                      if s.get("mem_total_gb") else self.t("res.na"))
        try:
            dp = estimator.GB * (s.get("disk_total_gb") or 0)
            dfree = estimator.GB * (s.get("disk_free_gb") or 0)
            pct = (100.0 * (dp - dfree) / dp) if dp else None
            self._set_res("res.disk", pct,
                          self.t("res.disk_val", free=fmt_bytes(dfree),
                                 total=fmt_bytes(dp)) if dp else self.t("res.na"))
        except Exception:
            self._set_res("res.disk", None, self.t("res.na"))
        rss = s.get("proc_rss_gb")
        self._set_res("res.rss", None,
                      ("%.1f GB" % rss) if rss else self.t("res.na"))
        self.after(1500, self._sample_resources)

    def _set_res(self, key, pct, text):
        bar, lab = self.res[key]
        if pct is None:
            bar.configure(value=0)
        else:
            bar.configure(value=max(0.0, min(100.0, float(pct))))
        lab.configure(text=text)

    # ------------------------------------------------------------ 轮询
    def _drain_ui(self):
        """后台线程的结果搬回主线程（Tk 控件只能在主线程改）。"""
        while True:
            try:
                fn, arg = self._ui_q.get_nowait()
            except queue.Empty:
                return
            try:
                fn(arg)
            except Exception:
                pass

    def _poll(self):
        self._drain_ui()
        if self._events_path and os.path.isfile(self._events_path):
            try:
                size = os.path.getsize(self._events_path)
            except OSError:
                size = 0
            if size < self._ev_off:      # 被重写（新一轮/外部重启）
                self.state = runlog.RunState()
                self._ev_off = 0
                self._log_data = []
                self._render_log(force=True)
                self._append_log(self.t("msg.resync"))
            recs, self._ev_off = runlog.read_events(self._events_path, self._ev_off)
            if recs:
                for r in recs:
                    self.state.feed(r)
                    if r.get("ev") == "log":
                        self._append_log(r.get("line") or "", r.get("stage"))
                    elif r.get("ev") == "msg":
                        self._append_log(self.t(r.get("key") or "",
                                                **(r.get("args") or {})))
                self._apply(recs)
                last = recs[-1]
                if last.get("ev") == "run.end":
                    self._on_run_end(last)
                if last.get("stage") in ("C", "P1", "P6b"):
                    self._refresh_db()
        if self._proc is not None and self._proc.poll() is not None:
            self._on_proc_exit()
        self.after(400, self._poll)

    def _on_run_end(self, rec):
        status = rec.get("status") or "ok"
        try:
            self.bell()
        except tk.TclError:
            pass
        self._flash_taskbar()
        if status in ("fail", "check-fail"):
            self._taskbar.set_state("error")
        elif status == "stopped":
            self._taskbar.set_state("paused")
        else:
            self._taskbar.set_progress(1.0, state="normal")
            self.after(2500, self._taskbar.clear)
        self._set_title_progress(None)
        self._set_status(self._status_text(status))
        if status in ("fail", "check-fail"):
            self._set_status(self._status_text(status), color=C["fail"])
        elif status == "ok":
            self._set_status(self._status_text(status), color=C["ok"])
        self.btn_start.configure(state="normal")
        self.btn_resume.configure(state="normal")
        self.btn_stop.configure(state="disabled")
        self._refresh_db()

    def _on_proc_exit(self):
        rc = self._proc.poll()
        # 进程退出前可能刚写了最后几行事件：补读一次
        if self._events_path and os.path.isfile(self._events_path):
            try:
                size = os.path.getsize(self._events_path)
            except OSError:
                size = 0
            if size < self._ev_off:      # 被重写（新一轮/外部重启）
                self.state = runlog.RunState()
                self._ev_off = 0
                self._log_data = []
                self._render_log(force=True)
                self._append_log(self.t("msg.resync"))
            recs, self._ev_off = runlog.read_events(self._events_path, self._ev_off)
            for r in recs:
                self.state.feed(r)
                if r.get("ev") == "log":
                    self._append_log(r.get("line") or "")
            if recs:
                self._apply([])
        self._proc = None
        if self._stdout_fh:
            try:
                self._stdout_fh.close()
            except OSError:
                pass
            self._stdout_fh = None
        if self.state.ended is None:
            status = "stopped" if rc == 130 else ("ok" if rc == 0 else "fail")
            self.state.feed({"ev": "run.end", "status": status,
                             "t": datetime.datetime.now().timestamp()})
            self._apply([])
            if rc not in (0, 130):
                self._dump_stdout_tail()
        self.btn_start.configure(state="normal")
        self.btn_resume.configure(state="normal")
        self.btn_stop.configure(state="disabled")
        if self._stopping:
            self._stopping = False
            self._set_status(self.t("msg.stopped", t=self.v_elapsed.get()))

    def _dump_stdout_tail(self):
        wd = self.v_workdir.get().strip()
        logs = os.path.join(wd, "00_logs")
        if not os.path.isdir(logs):
            return
        files = sorted(f for f in os.listdir(logs)
                       if f.startswith("gui_stdout_") and f.endswith(".log"))
        if not files:
            return
        try:
            with open(os.path.join(logs, files[-1]), encoding="utf-8",
                      errors="replace") as fh:
                tail = fh.readlines()[-40:]
        except OSError:
            return
        self._append_log("---- child stdout (tail) ----")
        for line in tail:
            self._append_log(line.rstrip("\n"))

    def _set_status(self, text, color=None):
        self.lb_status.configure(text=text or "", fg=color or C["muted"])

    # ------------------------------------------------------------ 启动/停止
    def _default_workdir(self):
        base = os.path.dirname(os.path.dirname(os.path.dirname(
            os.path.abspath(__file__))))          # ORFscape/
        stamp = datetime.datetime.now().strftime("%Y%m%d%H%M%S")
        return os.path.join(base, "orfscape_projects", "proj_" + stamp)

    def _start(self, mode):
        form = self._form()
        if not form["raws"]:
            messagebox.showwarning("ORFscape", self.t("msg.need_raw"))
            return
        if mode != "resume":       # 续跑时可能只想从中间某步接着跑
            n_raw, _bytes = model.raw_stats(form["raws"])
            if n_raw == 0:
                messagebox.showwarning("ORFscape", self.t("msg.no_raw_found"))
                return
        wd = self.v_workdir.get().strip()
        if mode == "resume" and wd:
            sid = model.resume_stage(self.state)
            if sid is None:
                messagebox.showinfo("ORFscape", self.t("msg.finished",
                                                     t=self.v_elapsed.get()))
                return
        else:
            sid = None
        if not wd:
            wd = self._default_workdir()
            self.v_workdir.set(wd)
        try:
            os.makedirs(wd, exist_ok=True)
        except OSError as exc:
            messagebox.showerror("ORFscape", "%s\n%s" % (self.t("msg.bad_dir"), exc))
            return

        args = model.build_run_args(form, wd, self.lang)
        if sid:
            args += ["--from", sid]
            if self.state.failed() is not None:
                # 失败阶段的产物标记可能已存在，不带 --force 会被直接跳过
                args += ["--force"]
        if self.v_confirm.get() and not self._confirm_plan(wd, form, sid):
            return
        pkg_root = os.path.dirname(os.path.dirname(os.path.dirname(
            os.path.abspath(__file__))))
        cmd = [sys.executable, "-m", "orfscape"] + args

        logs = os.path.join(wd, "00_logs")
        os.makedirs(logs, exist_ok=True)
        stamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        self._events_path = os.path.join(logs, "progress.jsonl")
        self._ev_off = 0
        self.state = runlog.RunState()
        self._log_data = []
        self._render_log(force=True)
        try:
            self._stdout_fh = open(os.path.join(logs, "gui_stdout_%s.log" % stamp),
                                   "w", encoding="utf-8", newline="\n")
        except OSError:
            self._stdout_fh = None

        env = dict(os.environ)
        env["PYTHONUNBUFFERED"] = "1"
        env["PYTHONIOENCODING"] = "utf-8"
        try:
            self._proc = subprocess.Popen(
                cmd, cwd=pkg_root, env=env,
                stdout=self._stdout_fh or subprocess.DEVNULL,
                stderr=subprocess.STDOUT)
        except OSError as exc:
            messagebox.showerror("ORFscape", str(exc))
            self._proc = None
            return
        self._taskbar.set_state("normal")
        self._taskbar.set_progress(0.0)
        self._run_t0 = datetime.datetime.now()
        self.btn_start.configure(state="disabled")
        self.btn_resume.configure(state="disabled")
        self.btn_stop.configure(state="normal")
        self._set_status(self.t("msg.started",
                                t=self._run_t0.strftime("%H:%M:%S")))
        self._save_cfg()

    # ------------------------------------------------------------ 运行计划
    def _plan_rows(self, workdir, form, start_sid=None):
        """本地生成运行计划（不启动进程）：返回 [(sid, 标题, 是否跳过, cmd)]。"""
        from .. import orchestrator
        ns = types.SimpleNamespace(
            raw=model.common_raw_root(form.get("raws") or []) or None,
            thread=form.get("threads"), workers=form.get("workers"),
            dbs=form.get("dbs"), min_cos="auto", reuse=None, cards=None,
            no_esm=False, species_name=form.get("species") or None,
            species_code=None, start="P0", end="P8c", skip="")
        ids = [s[0] for s in orchestrator.STEPS]
        lo = ids.index(start_sid) if start_sid in ids else 0
        rows = []
        for sid, sub, _chk, marker, _w in orchestrator.STEPS[lo:]:
            try:
                done = orchestrator.marker_exists(workdir, marker)
            except (ValueError, OSError):
                done = False
            try:
                cmd = " ".join(orchestrator.build_cmd(sid, sub, ns, workdir))
            except Exception:
                cmd = "(command unavailable)"
            rows.append((sid, orchestrator.TITLES.get(sid, sid), done, cmd))
        return rows

    def _build_plan_window(self, rows, est=None):
        """把计划渲染成窗口，返回 (窗口, 选择结果)。"""
        win = tk.Toplevel(self)
        win.title(self.t("plan.title"))
        win.geometry("820x460+%d+%d" % (self.winfo_rootx() + 80,
                                        self.winfo_rooty() + 70))
        win.configure(bg=C["card"])
        win.transient(self)

        head = tk.Frame(win, bg=C["card"])
        head.pack(fill="x", padx=14, pady=(12, 4))
        ttl = tk.Label(head, text="", bg=C["card"], fg=C["text"], font=FONT_H2)
        ttl.pack(side="left")
        self._bind_text(ttl, "text", "plan.title")
        est = est or {}
        if est and not est.get("error"):
            d = est["disk"]
            summary = self.t("plan.est", t=fmt_sec(est["total_sec"]),
                             d=fmt_bytes(d["total"]),
                             r="%.1f GB" % est["ram"]["need_gb"])
        else:
            summary = ""
        tk.Label(head, text=summary, bg=C["card"], fg=C["muted"],
                 font=FONT).pack(side="right")

        hint = tk.Label(win, text="", bg=C["card"], fg=C["hint"],
                        font=("Segoe UI", 8), anchor="w")
        hint.pack(fill="x", padx=14)
        self._bind_text(hint, "text", "plan.hint")

        wrap = tk.Frame(win, bg=C["border"])
        wrap.pack(fill="both", expand=True, padx=14, pady=(6, 8))
        txt = tk.Text(wrap, font=FONT_MONO, bd=0, wrap="char", bg=C["field"],
                      fg=C["text"], highlightthickness=0, padx=6, pady=4)
        txt.pack(side="left", fill="both", expand=True, padx=1, pady=1)
        sb = ttk.Scrollbar(wrap, orient="vertical", command=txt.yview)
        sb.pack(side="right", fill="y")
        txt.configure(yscrollcommand=sb.set)
        for sid, title, done, cmd in rows:
            tag = "skip" if done else "run"
            txt.insert("end", "[%s] %-4s %s\n"
                       % (self.t("plan.skip") if done else self.t("plan.run"),
                          sid, title), tag)
            txt.insert("end", "      " + cmd + "\n\n", ("cmd", tag))
        txt.tag_configure("skip", foreground=C["idle"])
        txt.tag_configure("run", foreground=C["text"])
        txt.tag_configure("cmd", foreground=C["muted"], font=("Consolas", 8))
        txt.configure(state="disabled")

        foot = tk.Frame(win, bg=C["card"])
        foot.pack(fill="x", padx=14, pady=(0, 12))
        choice = {"ok": False}

        def _ok(_ev=None):
            choice["ok"] = True
            win.destroy()

        b_ok = ttk.Button(foot, text="", width=12, command=_ok)
        b_ok.pack(side="right", padx=(8, 0))
        self._bind_text(b_ok, "text", "btn.start")
        b_no = ttk.Button(foot, text="", width=10, command=win.destroy)
        b_no.pack(side="right")
        self._bind_text(b_no, "text", "btn.cancel")
        win.bind("<Return>", _ok)
        win.bind("<Escape>", lambda e: win.destroy())
        return win, choice

    def _confirm_plan(self, workdir, form, start_sid=None):
        """开跑前把要执行的命令列出来；用户取消返回 False。"""
        try:
            rows = self._plan_rows(workdir, form, start_sid)
        except Exception:
            return True
        est = self._est
        if not est or est.get("error"):
            try:
                est = estimator.estimate(workdir, form.get("raws"),
                                         threads=form.get("threads") or 16,
                                         workers=form.get("workers") or 4,
                                         dbs=form.get("dbs") or "aug,6f")
            except Exception:
                est = None
        win, choice = self._build_plan_window(rows, est)
        win.grab_set()
        self.wait_window(win)
        return choice["ok"]

    def _stop(self):
        wd = self.v_workdir.get().strip()
        if not wd:
            return
        flag = os.path.join(wd, "00_logs", "stop.flag")
        try:
            os.makedirs(os.path.dirname(flag), exist_ok=True)
            with open(flag, "w", encoding="utf-8") as fh:
                fh.write("stop\n")
        except OSError:
            pass
        self._stopping = True
        self._set_status(self.t("msg.stop_req"), color=C["warn"])
        self.after(15000, self._force_kill)

    def _force_kill(self):
        if self._proc is not None and self._proc.poll() is None:
            try:
                if os.name == "nt":
                    subprocess.run(["taskkill", "/F", "/T", "/PID",
                                    str(self._proc.pid)], capture_output=True)
                else:
                    self._proc.terminate()
            except OSError:
                pass

    def _results_dir(self, wd):
        """结果目录：优先 P13_results，其次 08_report，最后项目根。"""
        for sub in ("P13_results", "08_report"):
            cand = os.path.join(wd, sub)
            if os.path.isdir(cand):
                return cand
        return wd

    def _open_results(self):
        wd = self.v_workdir.get().strip()
        if not wd or not os.path.isdir(wd):
            messagebox.showwarning("ORFscape", self.t("msg.need_dir"))
            return
        target = self._results_dir(wd)
        try:
            if os.name == "nt":
                os.startfile(target)                 # noqa: S606
            elif sys.platform == "darwin":
                subprocess.Popen(["open", target])
            else:
                subprocess.Popen(["xdg-open", target])
        except OSError as exc:
            messagebox.showerror("ORFscape", str(exc))

    def _export_summary(self):
        wd = self.v_workdir.get().strip()
        if not wd or not os.path.isdir(wd):
            messagebox.showwarning("ORFscape", self.t("msg.need_dir"))
            return
        if not self.state.order or not self.state.n_events:
            messagebox.showinfo("ORFscape", self.t("msg.no_log"))
            return
        try:
            path = model.export_summary(self.state, wd, est=self._est)
        except (OSError, ValueError) as exc:
            messagebox.showerror("ORFscape", str(exc))
            return
        self._set_status(self.t("msg.exported", path=os.path.basename(str(path))),
                         color=C["ok"])

    def _on_close(self):
        try:
            self._taskbar.close()
        except Exception:
            pass
        try:
            self.cfg["geometry"] = self.geometry()
        except tk.TclError:
            pass
        if self._proc is not None and self._proc.poll() is None:
            if not messagebox.askyesno(
                    "ORFscape", "A run is in progress. Stop it and quit?"):
                return
            self._stop()
        self._save_cfg()
        self.destroy()


def main(prefill_workdir=None):
    app = App(prefill_workdir=prefill_workdir)
    app.mainloop()
    return 0


if __name__ == "__main__":
    sys.exit(main())
