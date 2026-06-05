# Quickstart

## Reading a mesh

```python
import meshio

mesh = meshio.read("mesh.msh")
# or explicitly specify the format:
mesh = meshio.read("mesh.msh", file_format="gmsh")
```

`read` accepts a file path (string or `os.PathLike`) or an open file buffer. When a buffer is used, `file_format` is required.

After reading:

```python
mesh.points        # numpy array, shape (num_points, dim)
mesh.cells         # list of CellBlock objects
mesh.point_data    # dict of str -> numpy array
mesh.cell_data     # dict of str -> list of numpy arrays (one per CellBlock)
mesh.point_sets    # dict of str -> numpy array of point indices
mesh.cell_sets     # dict of str -> list of numpy arrays of cell indices
mesh.field_data    # dict of str -> numpy array (scalar metadata)
```

## Writing a mesh

```python
mesh.write("out.vtk")
# or with explicit format:
mesh.write("out.vtk", file_format="vtk")
```

## Constructing and writing from scratch

```python
import numpy as np
import meshio

points = np.array([
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [1.0, 1.0, 0.0],
    [0.0, 1.0, 0.0],
    [2.0, 0.0, 0.0],
    [2.0, 1.0, 0.0],
])
cells = [
    ("triangle", np.array([[0, 1, 2], [0, 2, 3]])),
    ("quad",     np.array([[1, 4, 5, 2]])),
]

mesh = meshio.Mesh(
    points,
    cells,
    point_data={"temperature": np.array([0.3, -1.2, 0.5, 0.7, 0.0, -3.0])},
    cell_data={"material": [np.array([1, 1]), np.array([2])]},
)
mesh.write("result.vtu")
```

The shorthand `write_points_cells` skips constructing a `Mesh` object:

```python
meshio.write_points_cells("result.vtu", points, cells,
                           point_data={"temperature": ...})
```

## Quick format conversion

```python
mesh = meshio.read("input.msh")
mesh.write("output.vtu")
```

Or from the command line:

```sh
meshio convert input.msh output.vtu
```

## Inspecting a mesh

```python
print(mesh)
# <meshio mesh object>
#   Number of points: 6
#   Number of cells:
#     triangle: 2
#     quad: 1
#   Point data: temperature
#   Cell data: material
```
