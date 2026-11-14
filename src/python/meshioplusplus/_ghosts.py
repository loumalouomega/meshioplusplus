"""The cell/point/field data names for ghost cells, owning parts and global ids.

The Python twin of ``src/cpp/src/detail/ghost_names.hpp``; see ``doc/ghosts.md``
for the convention. ``tests/python/test_ghost_names.py`` asserts the two stay
in step.
"""

#: cell_data, int: 0 = owned by its piece, L >= 1 = ghost at BFS layer depth L.
PARTITION_GHOST = "partition:ghost"

#: cell_data / point_data, int: owning part (0-based) of a ghost cell
#: (DUPLICATEPOINT owner for points).
PARTITION_OWNER = "partition:owner"

#: cell_data, int: the piece a cell belongs to (pvtu/pvtp/pvd/vtx read; merge output).
PARTITION_PART = "partition:part"

#: point_data / cell_data, int: 0-based input indices (partition, split).
PARTITION_ORIGINAL_POINT_ID = "partition:original_point_id"
PARTITION_ORIGINAL_CELL_ID = "partition:original_cell_id"

#: point_data / cell_data, int: vendor/global ids (radioss, dpf, refine, cgns,
#: exodus maps when non-identity).
GLOBAL_NODE_ID = "global:node_id"
GLOBAL_CELL_ID = "global:cell_id"

#: point_data / cell_data, uint8: VTK-native bitfield (DUPLICATECELL = 1 cells,
#: DUPLICATEPOINT = 1 points); kept verbatim by vtu/vtp/vtx/pvtu.
VTK_GHOST_TYPE = "vtkGhostType"

#: VTK-native global ids carried by vtu/vtp/vtx/pvtu.
VTK_ORIGINAL_POINT_IDS = "vtkOriginalPointIds"
VTK_ORIGINAL_CELL_IDS = "vtkOriginalCellIds"

#: cell_data, int: Gmsh-native 1/0 flag from ``$GhostElements``.
GMSH_GHOST = "gmsh:ghost"

#: VTK bit values of ``vtkGhostType``.
VTK_DUPLICATE_CELL = 1
VTK_DUPLICATE_POINT = 1
