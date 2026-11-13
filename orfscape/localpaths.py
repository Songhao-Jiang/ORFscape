# -*- coding: utf-8 -*-
"""本机外部工具路径解析：环境变量 > 本机配置 > 常见安装目录。

第三方软件（pFind、isoPSV、BLAST+、Miniconda…）位置因机器而异，不能写死在
源码里。解析顺序：

1. 环境变量 ``ORFSCAPE_<KEY>``（如 ``ORFSCAPE_ISOPSV``）；
2. 包目录下的 ``paths.local.json``（**不入库**，由 ``paths.local.example.json``
   复制后填写；``.gitignore`` 已排除）；
3. 调用方给的常见安装目录（存在才用）。

找不到时返回空串，调用方负责给出"先装/先配"的可读提示，绝不静默用错路径。
"""
import json
import os

_CACHE = None


def _load():
    """读一次 paths.local.json；损坏就当空表（不让配置问题冒泡成崩溃）。"""
    global _CACHE
    if _CACHE is not None:
        return _CACHE
    data = {}
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "paths.local.json")
    try:
        if os.path.isfile(path):
            with open(path, encoding="utf-8") as fh:
                loaded = json.load(fh)
            if isinstance(loaded, dict):
                data = loaded
    except (OSError, ValueError):
        data = {}
    _CACHE = data
    return data


def get(key, *fallbacks, **kw):
    """解析一个外部路径。

    key        配置键（``isopsv`` / ``pfind`` / ``blast`` / ``blastn`` …）
    fallbacks  候选目录，按顺序取第一个存在的
    env        覆盖环境变量名（默认 ``ORFSCAPE_<KEY 大写>``）
    """
    env = kw.get("env") or ("ORFSCAPE_" + key.upper())
    val = (os.environ.get(env) or "").strip()
    if val:
        return val
    val = str(_load().get(key, "") or "").strip()
    if val:
        return val
    for cand in fallbacks:
        if cand and os.path.isdir(cand):
            return cand
    return ""


def exe(key, *fallbacks, **kw):
    """同 :func:`get`，但回落到 PATH 里的可执行文件（``shutil.which``）。"""
    val = get(key, *fallbacks, **kw)
    if val:
        return val
    import shutil
    return shutil.which(kw.get("name") or key) or ""
