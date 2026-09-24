import os
from contextlib import contextmanager


def is_buffer(obj, mode):
    return ("r" in mode and hasattr(obj, "read")) or (
        "w" in mode and hasattr(obj, "write")
    )


@contextmanager
def open_file(path_or_buf, mode="r", **kwargs):
    if is_buffer(path_or_buf, mode):
        yield path_or_buf
    else:
        with open(path_or_buf, mode, **kwargs) as f:
            yield f


# Z88's input and output files have fixed names; the ``.txt`` extension alone
# would pick the xyz reader.
_Z88_FILENAMES = ("z88i1.txt", "z88structure.txt", "z88o2.txt", "z88o3.txt")


def is_z88_filename(path) -> bool:
    """Whether ``path``'s basename is one of Z88's fixed file names."""
    try:
        name = os.path.basename(os.fspath(path))
    except TypeError:
        return False
    return name.lower() in _Z88_FILENAMES
