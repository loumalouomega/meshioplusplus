"""A pytest plugin that makes ``import meshio`` import meshio++ instead.

Upstream meshio's own test suite is written against ``meshio``; with this plugin
loaded (``-p alias_plugin``, with this directory on ``PYTHONPATH``) every
``import meshio`` and ``import meshio.<sub>`` resolves to ``meshioplusplus`` and
the same submodule, so the suite exercises meshio++ through the drop-in name
(roadmap 3.18). Nothing is installed under the name ``meshio``.
"""

import importlib
import importlib.abc
import importlib.util
import sys

_UPSTREAM = "meshio"
_OURS = "meshioplusplus"


class _AliasLoader(importlib.abc.Loader):
    def __init__(self, target):
        self._target = target

    def create_module(self, spec):
        return importlib.import_module(self._target)

    def exec_module(self, module):  # already executed under its real name
        pass


class _AliasFinder(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname != _UPSTREAM and not fullname.startswith(_UPSTREAM + "."):
            return None
        real = _OURS + fullname[len(_UPSTREAM) :]
        try:
            importlib.import_module(real)
        except ImportError:
            return None
        return importlib.util.spec_from_loader(fullname, _AliasLoader(real))


def pytest_configure(config):
    if not any(isinstance(f, _AliasFinder) for f in sys.meta_path):
        sys.meta_path.insert(0, _AliasFinder())
