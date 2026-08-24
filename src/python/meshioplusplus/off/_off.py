"""
I/O for the OFF surface format, cf.
<https://en.wikipedia.org/wiki/OFF_(file_format)>,
<http://www.geomview.org/docs/html/OFF.html>.
"""

import numpy as np

from .._common import warn
from .._exceptions import ReadError
from .._files import open_file
from .._mesh import CellBlock, Mesh
from .._provenance import TAG as _PROVENANCE_TAG


def read(filename):
    with open_file(filename) as f:
        points, cells = read_buffer(f)
    return Mesh(points, cells)


def read_buffer(f):
    # assert that the first line reads `OFF`
    line = f.readline()

    if isinstance(line, (bytes, bytearray)):
        raise ReadError("Expected text buffer, not bytes.")

    if line.strip() != "OFF":
        raise ReadError("Expected the first line to be `OFF`.")

    # fast forward to the next significant line
    while True:
        line = f.readline().strip()
        if line and line[0] != "#":
            break

    # This next line contains:
    # <number of vertices> <number of faces> <number of edges>
    num_verts, num_faces, _ = line.split(" ")
    num_verts = int(num_verts)
    num_faces = int(num_faces)

    verts = np.fromfile(f, dtype=float, count=3 * num_verts, sep=" ").reshape(
        num_verts, 3
    )

    # Faces are grouped by vertex count into triangle (3), quad (4), or
    # polygon (else) blocks; a run of same-count faces stays in one block
    # until the count changes.
    cells = []
    run_n = None
    run_rows = []

    def flush():
        if not run_rows:
            return
        name = {3: "triangle", 4: "quad"}.get(run_n, "polygon")
        cells.append(CellBlock(name, np.array(run_rows, dtype=int)))

    for _ in range(num_faces):
        n = int(np.fromfile(f, dtype=int, count=1, sep=" ")[0])
        if n < 3:
            raise ReadError("OFF: faces must have at least 3 vertices")
        row = np.fromfile(f, dtype=int, count=n, sep=" ")
        if n != run_n:
            flush()
            run_n = n
            run_rows = []
        run_rows.append(row)
    flush()

    return verts, cells


def write(filename, mesh):
    if mesh.points.shape[1] == 2:
        warn(
            "OFF requires 3D points, but 2D points given. "
            "Appending 0 as third component."
        )
        points = np.column_stack([mesh.points, np.zeros_like(mesh.points[:, 0])])
    else:
        points = mesh.points

    face_blocks = [c for c in mesh.cells if c.type in ("triangle", "quad", "polygon")]
    skip = [c for c in mesh.cells if c.type not in ("triangle", "quad", "polygon")]
    if skip:
        string = ", ".join(item.type for item in skip)
        warn(f"OFF only supports triangle/quad/polygon cells. Skipping {string}.")

    num_faces = sum(len(c.data) for c in face_blocks)

    with open(filename, "wb") as fh:
        fh.write(b"OFF\n")
        fh.write(f"# {_PROVENANCE_TAG}\n\n".encode())

        # counts
        c = f"{mesh.points.shape[0]} {num_faces} {0}\n\n"
        fh.write(c.encode())

        # vertices
        # np.savetxt(fh, mesh.points, "%r")  # slower
        fmt = " ".join(["{}"] * points.shape[1])
        out = "\n".join([fmt.format(*row) for row in points]) + "\n"
        fh.write(out.encode())

        # faces (each block may be a uniform ndarray or, for a ragged
        # "polygon" block, a Python list of per-face node arrays)
        lines = []
        for block in face_blocks:
            for row in block.data:
                lines.append(f"{len(row)} " + " ".join(str(i) for i in row))
        if lines:
            fh.write(("\n".join(lines) + "\n").encode())


# NOTE: format registration now lives in meshioplusplus/off/__init__.py, which wraps the
# reader/writer below with the C++-backed fast paths.
