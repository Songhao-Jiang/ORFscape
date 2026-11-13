# -*- coding: utf-8 -*-
"""配色与主题。

- 所有颜色都走 token（`C["card"]`），不在控件里写字面量；
- 切换主题只做两件事：`C` 原地更新 + 按"旧色 → 新色"映射重涂控件；
- 映射表由两套调色板自动生成，不需要为每个控件登记角色。

    from .theme import C, set_theme, color_map
    set_theme("dark")
"""

LIGHT = {
    "bg": "#F5F6F8",
    "card": "#FFFFFF",
    "border": "#DCE0E6",
    "text": "#1F2328",
    "muted": "#6B7280",
    "hint": "#9AA1AA",
    "accent": "#2563EB",
    "ok": "#16A34A",
    "warn": "#D97706",
    "fail": "#DC2626",
    "idle": "#9CA3AF",
    "chip": "#F3F4F6",
    "chip_run": "#E8F0FE",
    "chip_ok": "#EAF7EE",
    "chip_bad": "#FDECEC",
    "track": "#EDEFF2",
    "field": "#FFFFFF",
    "zebra": "#FAFBFC",
    "sel": "#E8F0FE",
    "hit": "#FFE8A3",
    "find_bg": "#FFD7A8",
    "log_bg": "#FBFBFC",
    "log_fg": "#374151",
    "ok_border": "#BFE3CB",
    "bad_border": "#F3C0C0",
    "stop_bg": "#FDF6EC",
    "stop_border": "#F0D9B5",
}

DARK = {
    "bg": "#16181D",
    "card": "#1E2128",
    "border": "#30343D",
    "text": "#E8EAED",
    "muted": "#A0A8B4",
    "hint": "#7A828E",
    "accent": "#5B9DFF",
    "ok": "#46C46A",
    "warn": "#E3A33B",
    "fail": "#FF6B64",
    "idle": "#6C7480",
    "chip": "#24272E",
    "chip_run": "#1B2A41",
    "chip_ok": "#1B2C22",
    "chip_bad": "#35211F",
    "track": "#2C3038",
    "field": "#131519",
    "zebra": "#21242B",
    "sel": "#26354D",
    "hit": "#4A3A00",
    "find_bg": "#5A4300",
    "log_bg": "#121418",
    "log_fg": "#C9D1D9",
    "ok_border": "#2E5338",
    "bad_border": "#5E2B27",
    "stop_bg": "#3A2E1B",
    "stop_border": "#6E4E20",
}

THEMES = {"light": LIGHT, "dark": DARK}
C = dict(LIGHT)          # 当前生效的调色板（原地更新，全局共享）
_current = "light"


def get_theme():
    return _current


def set_theme(name):
    """切换调色板；返回实际生效的主题名。"""
    global _current
    name = str(name or "").lower()
    if name not in THEMES:
        name = "light"
    _current = name
    C.clear()
    C.update(THEMES[name])
    return name


def color_map():
    """生成 "旧色 → 新色" 映射（当前主题 ↔ 另一套）。

    同一色值被多个 token 共用时（如 card 与 field 都是白），按调色板里的
    书写顺序先到先得；控件角色上的特例（Entry 用 field）由调用方单独覆盖。
    """
    other = DARK if _current == "light" else LIGHT
    out = {}
    for key, val in C.items():          # dict 顺序即优先级
        tgt = other.get(key)
        if not tgt or tgt.lower() == val.lower():
            continue
        out.setdefault(val.lower(), tgt)     # 单向：旧 → 新
    return out


def map_color(value, mapping):
    """把控件当前的色值按映射换成另一套；不认识的原样返回。"""
    if not value:
        return value
    low = str(value).strip().lower()
    return mapping.get(low, value)
