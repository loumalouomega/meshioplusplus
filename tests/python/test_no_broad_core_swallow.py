"""Guards the native-fallback shims from regressing to a blanket ``except Exception``.

Every format package wraps its ``_core`` call in a shim that falls back to the
Python reference implementation. The bare ``except Exception: pass`` that used
to do it hid why the C++ path declined and swallowed real bugs. The decision now
lives in ``meshioplusplus._fallback.core_declined``; this static guard fails any
handler that catches ``Exception`` around a ``_core.`` call without asking it.

Scoped to the per-format packages (``<fmt>/*.py``). The package-root operation
shims (``_clean.py``, ``_data_average.py``, ...) have the same defect shape but
are a separate scope.
"""

import ast
import pathlib

REPO = pathlib.Path(__file__).resolve().parent.parent.parent
PKG = REPO / "src" / "python" / "meshioplusplus"

_SKIP_DIRS = {"_cli", "mcp", "_viewer_assets", "__pycache__"}


def _format_python_files():
    for path in sorted(PKG.glob("*/*.py")):
        if path.parent.name in _SKIP_DIRS:
            continue
        yield path


def _calls_core(node):
    return any(
        isinstance(n, ast.Call)
        and isinstance(n.func, ast.Attribute)
        and isinstance(n.func.value, ast.Name)
        and n.func.value.id == "_core"
        for n in ast.walk(node)
    )


def _catches_everything(handler):
    t = handler.type
    if t is None:
        return True
    names = t.elts if isinstance(t, ast.Tuple) else [t]
    return any(
        isinstance(n, ast.Name) and n.id in ("Exception", "BaseException")
        for n in names
    )


def _asks_core_declined(handler):
    return any(
        isinstance(n, ast.Call)
        and isinstance(n.func, ast.Name)
        and n.func.id == "core_declined"
        for n in ast.walk(handler)
    )


def test_no_broad_except_around_core_without_core_declined():
    violations = []
    for path in _format_python_files():
        tree = ast.parse(path.read_text())
        for node in ast.walk(tree):
            if not isinstance(node, ast.Try):
                continue
            if not any(_calls_core(stmt) for stmt in node.body):
                continue
            for handler in node.handlers:
                if _catches_everything(handler) and not _asks_core_declined(handler):
                    violations.append(f"{path.relative_to(REPO)}:{handler.lineno}")

    assert not violations, (
        "a broad `except` around a `_core.` call must ask "
        "`meshioplusplus._fallback.core_declined` whether to fall back, "
        "not swallow the reason:\n" + "\n".join(violations)
    )


def test_the_guard_actually_sees_the_shims():
    # A guard that matches nothing passes vacuously; pin that it finds the
    # 47 shim packages' handlers.
    seen = 0
    for path in _format_python_files():
        tree = ast.parse(path.read_text())
        for node in ast.walk(tree):
            if isinstance(node, ast.Try) and any(_calls_core(s) for s in node.body):
                seen += sum(1 for h in node.handlers if _catches_everything(h))
    assert seen >= 80, seen
