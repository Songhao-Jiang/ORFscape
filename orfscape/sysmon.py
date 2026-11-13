# -*- coding: utf-8 -*-
"""系统资源采样（纯标准库，Windows 优先，Linux 降级）。

GUI 需要的是"这台机器吃不吃得下"：总 CPU 占用、内存占用、目标盘剩余，
以及当前流水线进程**整棵进程树**的常驻内存（python 自己 + pFind 子进程）。
psutil 不一定装了，所以这里用 ctypes / /proc 直接问内核；任何一项取不到都
返回 None，界面显示 n/a，不影响主流程。

    s = Sampler()
    s.sample(pid)   # {"cpu_pct":63.2, "mem_used_gb":.., "proc_rss_gb":.., "proc_tree_n":7}
"""
import os
import shutil
import sys

_WIN = sys.platform.startswith("win")
_LIN = sys.platform.startswith("linux")


# ---------------------------------------------------------------- Windows
if _WIN:
    import ctypes
    from ctypes import wintypes

    class _FILETIME(ctypes.Structure):
        _fields_ = [("dwLowDateTime", wintypes.DWORD),
                    ("dwHighDateTime", wintypes.DWORD)]

    class _MEMSTATUS(ctypes.Structure):
        _fields_ = [("dwLength", wintypes.DWORD),
                    ("dwMemoryLoad", wintypes.DWORD),
                    ("ullTotalPhys", ctypes.c_ulonglong),
                    ("ullAvailPhys", ctypes.c_ulonglong),
                    ("ullTotalPageFile", ctypes.c_ulonglong),
                    ("ullAvailPageFile", ctypes.c_ulonglong),
                    ("ullTotalVirtual", ctypes.c_ulonglong),
                    ("ullAvailVirtual", ctypes.c_ulonglong),
                    ("ullAvailExtendedVirtual", ctypes.c_ulonglong)]

    class _PMC(ctypes.Structure):
        _fields_ = [("cb", wintypes.DWORD), ("PageFaultCount", wintypes.DWORD),
                    ("PeakWorkingSetSize", ctypes.c_size_t),
                    ("WorkingSetSize", ctypes.c_size_t),
                    ("QuotaPeakPagedPoolUsage", ctypes.c_size_t),
                    ("QuotaPagedPoolUsage", ctypes.c_size_t),
                    ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t),
                    ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
                    ("PagefileUsage", ctypes.c_size_t),
                    ("PeakPagefileUsage", ctypes.c_size_t)]

    class _PROCESSENTRY32(ctypes.Structure):
        _fields_ = [("dwSize", wintypes.DWORD),
                    ("cntUsage", wintypes.DWORD),
                    ("th32ProcessID", wintypes.DWORD),
                    ("th32DefaultHeapID", ctypes.c_size_t),
                    ("th32ModuleID", wintypes.DWORD),
                    ("cntThreads", wintypes.DWORD),
                    ("th32ParentProcessID", wintypes.DWORD),
                    ("pcPriClassBase", wintypes.LONG),
                    ("dwFlags", wintypes.DWORD),
                    ("szExeFile", ctypes.c_char * 260)]

    _k32 = ctypes.WinDLL("kernel32", use_last_error=True)
    _psapi = ctypes.WinDLL("psapi", use_last_error=True)
    _TH32CS_SNAPPROCESS = 0x00000002
    _INVALID_HANDLE = ctypes.c_void_p(-1).value

    def _ft(v):
        return (v.dwHighDateTime << 32) | v.dwLowDateTime

    def _cpu_times():
        idle, kern, user = _FILETIME(), _FILETIME(), _FILETIME()
        if not _k32.GetSystemTimes(ctypes.byref(idle), ctypes.byref(kern),
                                   ctypes.byref(user)):
            return None
        return _ft(idle), _ft(kern), _ft(user)

    def _mem():
        st = _MEMSTATUS()
        st.dwLength = ctypes.sizeof(_MEMSTATUS)
        if not _k32.GlobalMemoryStatusEx(ctypes.byref(st)):
            return None
        return (st.ullTotalPhys - st.ullAvailPhys, st.ullTotalPhys)

    def _proc_rss(pid):
        if not pid:
            return None
        PROCESS_QUERY_LIMITED_INFORMATION = 0x1000
        h = _k32.OpenProcess(PROCESS_QUERY_LIMITED_INFORMATION, False, int(pid))
        if not h:
            return None
        try:
            pmc = _PMC()
            pmc.cb = ctypes.sizeof(_PMC)
            if not _psapi.GetProcessMemoryInfo(h, ctypes.byref(pmc), pmc.cb):
                return None
            return float(pmc.WorkingSetSize)
        finally:
            _k32.CloseHandle(h)

    def _children_map():
        """{ppid: [pid, ...]}（Toolhelp32 快照）。"""
        snap = _k32.CreateToolhelp32Snapshot(_TH32CS_SNAPPROCESS, 0)
        if not snap or snap == _INVALID_HANDLE:
            return {}
        out = {}
        try:
            ent = _PROCESSENTRY32()
            ent.dwSize = ctypes.sizeof(_PROCESSENTRY32)
            ok = _k32.Process32First(snap, ctypes.byref(ent))
            while ok:
                out.setdefault(int(ent.th32ParentProcessID), []).append(
                    int(ent.th32ProcessID))
                ok = _k32.Process32Next(snap, ctypes.byref(ent))
        finally:
            _k32.CloseHandle(snap)
        return out

    def _proc_tree_rss(pid):
        """整棵进程树（含 pFind 等孙进程）的 RSS 之和。"""
        pid = int(pid)
        kids = _children_map()
        seen, total, n = set(), 0.0, 0
        stack = [pid]
        while stack:
            cur = stack.pop()
            if cur in seen:
                continue
            seen.add(cur)
            rss = _proc_rss(cur)
            if rss:
                total += rss
                n += 1
            stack.extend(kids.get(cur, ()))
        return (total, n) if n else (None, 0)

else:
    def _cpu_times():
        try:
            with open("/proc/stat") as fh:
                f = fh.readline().split()
            vals = [int(x) for x in f[1:]]
            idle = vals[3] + (vals[4] if len(vals) > 4 else 0)
            total = sum(vals[:8])
            return idle, total, 0
        except OSError:
            return None

    def _mem():
        try:
            info = {}
            with open("/proc/meminfo") as fh:
                for line in fh:
                    k, _, v = line.partition(":")
                    info[k.strip()] = int(v.split()[0]) * 1024
            total = info.get("MemTotal")
            avail = info.get("MemAvailable", info.get("MemFree"))
            if not total or not avail:
                return None
            return (total - avail, total)
        except OSError:
            return None

    def _proc_rss(pid):
        try:
            with open("/proc/%d/status" % int(pid)) as fh:
                for line in fh:
                    if line.startswith("VmRSS:"):
                        return float(line.split()[1]) * 1024
        except (OSError, TypeError, ValueError):
            pass
        return None

    def _children_map():
        out = {}
        try:
            names = os.listdir("/proc")
        except OSError:
            return out
        for name in names:
            if not name.isdigit():
                continue
            try:
                with open("/proc/%s/stat" % name) as fh:
                    data = fh.read()
            except OSError:
                continue
            # comm 里可能有空格/括号，从右括号之后切
            close = data.rfind(")")
            if close < 0:
                continue
            fields = data[close + 2:].split()
            if len(fields) < 2:
                continue
            try:
                out.setdefault(int(fields[1]), []).append(int(name))
            except ValueError:
                continue
        return out

    def _proc_tree_rss(pid):
        pid = int(pid)
        kids = _children_map()
        seen, total, n = set(), 0.0, 0
        stack = [pid]
        while stack:
            cur = stack.pop()
            if cur in seen:
                continue
            seen.add(cur)
            rss = _proc_rss(cur)
            if rss:
                total += rss
                n += 1
            stack.extend(kids.get(cur, ()))
        return (total, n) if n else (None, 0)

    _LIN = True


GB = 1024.0 ** 3


class Sampler(object):
    """差值式 CPU 采样 + 即时内存/磁盘/进程树 RSS。"""

    def __init__(self, disk_path=None):
        self.disk_path = disk_path or os.getcwd()
        self._prev = _cpu_times()

    def cpu_pct(self):
        cur = _cpu_times()
        if not cur or not self._prev:
            self._prev = cur
            return None
        pi, pk, pu = self._prev
        ci, ck, cu = cur
        self._prev = cur
        if _WIN:
            d_idle, d_tot = ci - pi, (ck - pk) + (cu - pu)
        else:
            d_idle, d_tot = ci - pi, ck - pk
        if d_tot <= 0:
            return None
        return max(0.0, min(100.0, 100.0 * (1.0 - float(d_idle) / d_tot)))

    def sample(self, pid=None):
        mem = _mem()
        rss, n_proc = (None, 0)
        if pid:
            try:
                rss, n_proc = _proc_tree_rss(pid)
            except Exception:
                rss, n_proc = None, 0
        out = {
            "cpu_pct": self.cpu_pct(),
            "cpu_n": os.cpu_count(),
            "mem_used_gb": round(mem[0] / GB, 2) if mem else None,
            "mem_total_gb": round(mem[1] / GB, 2) if mem else None,
            "mem_pct": round(100.0 * mem[0] / mem[1], 1) if mem and mem[1] else None,
            "proc_rss_gb": round(rss / GB, 2) if rss else None,
            "proc_tree_n": n_proc,
        }
        try:
            du = shutil.disk_usage(self.disk_path)
            out.update(disk_free_gb=round(du.free / GB, 2),
                       disk_total_gb=round(du.total / GB, 2),
                       disk_pct=round(100.0 * du.used / du.total, 1))
        except OSError:
            out.update(disk_free_gb=None, disk_total_gb=None, disk_pct=None)
        return out

    def set_disk_path(self, p):
        if p and os.path.isdir(p):
            self.disk_path = p


def snapshot(disk_path=None, pid=None):
    """一次性快照（无 CPU 差值，cpu_pct 为 None）。"""
    mem = _mem()
    out = {"cpu_pct": None, "cpu_n": os.cpu_count(),
           "mem_used_gb": round(mem[0] / GB, 2) if mem else None,
           "mem_total_gb": round(mem[1] / GB, 2) if mem else None,
           "mem_pct": round(100.0 * mem[0] / mem[1], 1) if mem and mem[1] else None,
           "proc_rss_gb": None, "proc_tree_n": 0}
    try:
        du = shutil.disk_usage(disk_path or os.getcwd())
        out.update(disk_free_gb=round(du.free / GB, 2),
                   disk_total_gb=round(du.total / GB, 2),
                   disk_pct=round(100.0 * du.used / du.total, 1))
    except OSError:
        out.update(disk_free_gb=None, disk_total_gb=None, disk_pct=None)
    return out
