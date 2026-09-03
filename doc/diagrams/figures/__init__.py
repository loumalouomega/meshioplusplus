"""The figure registry: name -> callable returning SVG text.

Order matters only for ``--list``; every entry is also embedded by at least
one page under ``doc/`` or by ``README.md`` (``gen_diagrams.py --check``
fails otherwise).
"""

from collections import OrderedDict

from . import architecture, cells, flows, refine, topology

REGISTRY = OrderedDict()
for _module in (architecture, flows, cells, topology, refine):
    REGISTRY.update(_module.FIGURES)
