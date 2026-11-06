from .._helpers import register_format
from .._parallel_xml import make_parallel_format
from ._pvtu import read as _py_read
from ._pvtu import write as _py_write
from ._pvtu import write_pieces as _py_write_pieces

read, write, write_pieces = make_parallel_format(
    "pvtu", "PUnstructuredGrid", ".vtu", _py_read, _py_write, _py_write_pieces
)

register_format("pvtu", [".pvtu"], read, {"pvtu": write})

__all__ = ["read", "write", "write_pieces"]
