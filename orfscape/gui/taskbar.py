# -*- coding: utf-8 -*-
"""Windows 任务栏进度（`ITaskbarList3`），纯 ctypes，失败一律降级为无操作。

用法::

    tb = TaskbarProgress()
    tb.bind(root)              # 传任意 tk 控件，取它的根窗口句柄
    tb.set_progress(0.42)      # 任务栏图标上出现 42%
    tb.set_state("error")      # 失败：红色
    tb.clear()                 # 清掉

约束：任何一步出错（非 Windows、COM 初始化失败、接口不受支持、句柄失效）
都置 ``self.ok = False`` 并静默返回，绝不抛到界面线程。
"""
import ctypes
import sys

_WIN = sys.platform.startswith("win")

# TBPF_*
STATES = {
    "none": 0x0,
    "indeterminate": 0x1,
    "normal": 0x2,
    "error": 0x4,
    "paused": 0x8,
}

GA_ROOT = 2
CLSCTX_INPROC_SERVER = 0x1

# vtable 槽位（IUnknown 0-2 + ITaskbarList 3-7 + ITaskbarList2 8 + ITaskbarList3 9-）
SLOT_RELEASE = 2
SLOT_HRINIT = 3
SLOT_SET_PROGRESS_VALUE = 9
SLOT_SET_PROGRESS_STATE = 10


class _GUID(ctypes.Structure):
    _fields_ = [("Data1", ctypes.c_ulong),
                ("Data2", ctypes.c_ushort),
                ("Data3", ctypes.c_ushort),
                ("Data4", ctypes.c_ubyte * 8)]


class TaskbarProgress(object):
    def __init__(self, enabled=True):
        self.ok = False
        self._ptr = None
        self._widget = None
        self._last = None          # (percent, state) 去抖
        if not (_WIN and enabled):
            return
        try:
            self._init_com()
        except Exception:
            self.ok = False

    # ---------------------------------------------------------------- 初始化
    def _init_com(self):
        ole32 = ctypes.windll.ole32
        ole32.CoInitializeEx.argtypes = [ctypes.c_void_p, ctypes.c_ulong]
        ole32.CoInitializeEx.restype = ctypes.c_long
        hr = ole32.CoInitializeEx(None, 0x2)      # COINIT_APARTMENTTHREADED
        if hr not in (0, 1, -2147417850):         # S_OK / S_FALSE / CHANGED_MODE
            return

        clsid, iid = _GUID(), _GUID()
        ole32.CLSIDFromString.argtypes = [ctypes.c_wchar_p,
                                          ctypes.POINTER(_GUID)]
        ole32.IIDFromString.argtypes = [ctypes.c_wchar_p,
                                        ctypes.POINTER(_GUID)]
        if ole32.CLSIDFromString("{56FDF344-FD6D-11D0-958A-006097C9A090}",
                                 ctypes.byref(clsid)) != 0:
            return
        if ole32.IIDFromString("{EA1AFB91-9E28-4B86-90E9-9E9F8A5EEA84}",
                               ctypes.byref(iid)) != 0:
            return

        ptr = ctypes.c_void_p()
        ole32.CoCreateInstance.argtypes = [ctypes.POINTER(_GUID),
                                           ctypes.c_void_p, ctypes.c_ulong,
                                           ctypes.POINTER(_GUID),
                                           ctypes.POINTER(ctypes.c_void_p)]
        hr = ole32.CoCreateInstance(ctypes.byref(clsid), None,
                                    CLSCTX_INPROC_SERVER, ctypes.byref(iid),
                                    ctypes.byref(ptr))
        if hr != 0 or not ptr:
            return
        self._ptr = ptr
        self._call(SLOT_HRINIT, ctypes.c_long)
        self.ok = True

    def _call(self, slot, restype, argtypes=(), args=()):
        """按 vtable 槽位调用；任何异常都降级为 None。"""
        if not self._ptr:
            return None
        try:
            vtbl = ctypes.cast(
                self._ptr, ctypes.POINTER(ctypes.POINTER(ctypes.c_void_p))
            ).contents
            proto = ctypes.WINFUNCTYPE(restype, ctypes.c_void_p, *argtypes)
            fn = proto(vtbl[slot])
            return fn(self._ptr, *args)
        except Exception:
            self.ok = False
            return None

    # ---------------------------------------------------------------- 对外
    def bind(self, widget):
        self._widget = widget

    def hwnd(self):
        if not self._widget:
            return 0
        try:
            return int(ctypes.windll.user32.GetAncestor(
                self._widget.winfo_id(), GA_ROOT))
        except Exception:
            return 0

    def set_state(self, name):
        if not self.ok:
            return
        if self._last and self._last[1] == name:
            return
        pct = self._last[0] if self._last else 0
        self._last = (pct, name)
        self._call(SLOT_SET_PROGRESS_STATE, ctypes.c_long,
                   (ctypes.c_void_p, ctypes.c_int),
                   (self.hwnd(), STATES.get(name, 0)))

    def set_progress(self, frac, state=None):
        """frac 0..1；百分比没变（且未指定状态）就不重复调用 COM。"""
        if not self.ok:
            return
        try:
            frac = max(0.0, min(1.0, float(frac)))
        except (TypeError, ValueError):
            return
        pct = int(round(frac * 100))
        prev_state = self._last[1] if self._last else "normal"
        if self._last and self._last[0] == pct and state is None:
            return
        self._last = (pct, state or prev_state)
        if state:
            self.set_state(state)
        self._call(SLOT_SET_PROGRESS_VALUE, ctypes.c_long,
                   (ctypes.c_void_p, ctypes.c_ulonglong, ctypes.c_ulonglong),
                   (self.hwnd(), pct, 100))

    def clear(self):
        if not self.ok:
            return
        self._last = None
        self._call(SLOT_SET_PROGRESS_STATE, ctypes.c_long,
                   (ctypes.c_void_p, ctypes.c_int), (self.hwnd(), 0))

    def close(self):
        try:
            self.clear()
        except Exception:
            pass
        if self._ptr:
            self._call(SLOT_RELEASE, ctypes.c_ulong, (), ())
            self._ptr = None
        self.ok = False
