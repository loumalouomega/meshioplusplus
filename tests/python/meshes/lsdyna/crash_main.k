*KEYWORD
$ crash_main.k: a small multi-file, mixed-format deck
*TITLE
crash test fixture
*INCLUDE_PATH
elements
*INCLUDE
crash_nodes.k
*INCLUDE
solids.k
*INCLUDE
shells.k
*INCLUDE
materials_not_shipped.k
*ELEMENT_BEAM
       1       5       1       2      27
       2       5       2       3      27
       3       5       4       5      27
*PART
matrix
         1         1         1
*PART
skin
         2         1         1
*PART
plate
         3         1         1
*PART
loose tets
         4         1         1
*PART
rods
         5         1         1
*SET_NODE_LIST_GENERATE
         1
         1         9         0         0
*SET_SHELL_LIST_TITLE
bottom skin
         2
         1         2         3         4
*SET_PART_LIST_TITLE
solid parts
         3
         1         3         4
*SET_SEGMENT_TITLE
loaded face
         4
         1         2         5         4
*CONTROL_TERMINATION
       1.0
*END
