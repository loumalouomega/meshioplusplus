from .._helpers import register_format
from .._parallel_xml import make_parallel_format
from ._pvtp import read as _py_read
from ._pvtp import write as _py_write
from ._pvtp import write_pieces as _py_write_pieces

read, write, write_pieces = make_parallel_format(
    "pvtp", "PPolyData", ".vtp", _py_read, _py_write, _py_write_pieces
)

register_format("pvtp", [".pvtp"], read, {"pvtp": write})

__all__ = ["read", "write", "write_pieces"]
