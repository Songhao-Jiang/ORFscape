# -*- coding: utf-8 -*-
"""轻量 lint：找未使用的 import（不装 pyflakes，纯标准库 AST）。

    python ORFscape/tests/test_gui_lint.py [--fix-hints]

只扫 GUI 与编排相关模块；有未使用 import 时退出码 1，并打印 `文件:行 名字`。
"""
import ast
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PKG = HERE.parent / "orfscape"
TARGETS = [
    PKG / "gui" / "app.py",
    PKG / "gui" / "model.py",
    PKG / "gui" / "i18n.py",
    PKG / "gui" / "theme.py",
    PKG / "gui" / "taskbar.py",
    PKG / "progress.py",
    PKG / "runlog.py",
    PKG / "estimator.py",
    PKG / "checks.py",
    PKG / "sysmon.py",
    PKG / "orchestrator.py",
    PKG / "core" / "eukbuild.py",
    PKG / "core" / "eukpresets.py",
    PKG / "core" / "netguard.py",
]


def bound_names(node):
    """一条 import 语句绑定到本模块的名字。"""
    out = []
    if isinstance(node, ast.Import):
        for a in node.names:
            out.append((a.asname or a.name.split(".")[0], node.lineno))
    elif isinstance(node, ast.ImportFrom):
        for a in node.names:
            if a.name == "*":
                continue
            out.append((a.asname or a.name, node.lineno))
    return out


def used_names(tree):
    used = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Name):
            used.add(node.id)
        elif isinstance(node, ast.Attribute):
            cur = node
            while isinstance(cur, ast.Attribute):
                cur = cur.value
            if isinstance(cur, ast.Name):
                used.add(cur.id)
    return used


def check(path):
    src = path.read_text(encoding="utf-8")
    tree = ast.parse(src)
    used = used_names(tree)
    problems = []
    for node in ast.walk(tree):
        if isinstance(node, (ast.Import, ast.ImportFrom)):
            for name, line in bound_names(node):
                if name not in used:
                    problems.append((path.name, line, name))
    return problems


def main():
    bad = []
    for path in TARGETS:
        if not path.is_file():
            continue
        bad.extend(check(path))
    if bad:
        for f, line, name in bad:
            print("unused import: %s:%d %s" % (f, line, name))
        print("%d unused import(s)" % len(bad))
        return 1
    print("OK — no unused imports in %d file(s)" % len(TARGETS))
    return 0


if __name__ == "__main__":
    sys.exit(main())
