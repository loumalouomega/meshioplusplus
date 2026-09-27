# Periodic node pairs

`match_periodic_nodes(mesh, slave, master, translate=…)` pairs every node of one boundary region with the node of another that a transform maps it onto — the pairs a periodic boundary condition ties together, whether that is a Kratos periodic condition, a Gmsh `$Periodic` section or a cyclic-symmetry constraint.

```python
import meshioplusplus as mp

mesh = mp.read("channel.msh")

# Opposite faces of a box, 2.0 apart in x:
pairs = mp.match_periodic_nodes(mesh, "inlet", "outlet", translate=(2.0, 0.0, 0.0))
pairs[:, 0]   # slave node ids, ascending
pairs[:, 1]   # the master node each one maps onto

# The cut faces of a 60-degree sector of a rotor:
pairs, report = mp.match_periodic_nodes(
    mesh, "cut_a", "cut_b", rotate=("z", 60.0), origin=(0, 0, 0), return_report=True
)
report["num_fixed"]      # nodes on the axis, which map onto themselves
report["max_residual"]   # how far the worst transformed node was from its master
```

## Which nodes

A region contributes the nodes it names: a point region its entries, a cell region the nodes of its cells, a side region the nodes of its facets — so a Gmsh physical surface, an Abaqus `*NSET` and an Exodus side set all work as they come. A region is named as for [`edit_regions`](/regions#editing-regions): by name, or by `{"name", "kind", "dim", "tag"}` when a name is shared.

## The transform

It maps a slave node's position onto its master's. It is given in [`transform`](/transform)'s vocabulary: `translate=(dx, dy, dz)`, `rotate=(axis, degrees)` with `axis` one of `"x"`, `"y"`, `"z"` or a 3-vector, turning about `origin` and applied before the translation, or `matrix=` a row-major 4x4 affine matrix (the layout Gmsh's `$Periodic` affine uses) overriding both.

## Guarantees

The nearest master node within `atol` wins, ties going to the lower node id. Two slave nodes landing on one master node is an error — the tolerance is too loose, or the regions overlap — and so is a slave node with no master, unless `require_complete=False`, in which case it is reported in `report["unmatched"]`. A slave node that the transform leaves in place and that is also a master node, such as a node on a rotation axis shared by both cut faces, is a *fixed point*: counted, not paired. The pairs are ordered by slave node id and are the same on every backend and thread count.

## Writing a Gmsh `$Periodic` section

The Python Gmsh writer emits `mesh.gmsh_periodic`, a list of `[dim, (slave_tag, master_tag), affine, node_pairs]` records ([Gmsh](/formats/gmsh)). One record from matched pairs:

```python
m = mp._periodic.periodic_matrix(translate=(2.0, 0.0, 0.0))
mesh.gmsh_periodic = [[2, (slave_tag, master_tag), m.reshape(-1), pairs]]
mp.write("channel_periodic.msh", mesh, file_format="gmsh")
```

## CLI

```sh
meshioplusplus periodic IN --slave NAME --master NAME \
    (--translate DX,DY,DZ | --rotate AXIS,DEG [--origin X,Y,Z] | --matrix M00,...,M33) \
    [--atol TOL] [--allow-incomplete] [--output pairs.csv] [--json]
```

Numbers are comma-separated; write a leading minus as `--translate=-2,0,0`. `--output` writes `slave,master` CSV rows of 0-based node ids.

## Other languages

```c
mio_region_selector s, t;
mio_region_selector_init(&s, "inlet");
mio_region_selector_init(&t, "outlet");
mio_periodic_opts opts;
mio_periodic_opts_init(&opts);          /* identity, atol 1e-8 */
opts.matrix[3] = 2.0;                   /* translate x */
mio_periodic_pairs* p = mio_match_periodic_nodes(mesh, &s, &t, &opts);
int64_t n;
const int64_t* slave = mio_periodic_pairs_slave(p, &n);
const int64_t* master = mio_periodic_pairs_master(p, NULL);
mio_periodic_pairs_free(p);
```

```fortran
call m%match_periodic_nodes('inlet', 'outlet', slave_ids, master_ids, translate=[2.0_real64, 0.0_real64, 0.0_real64])
```

```julia
p = match_periodic_nodes(mesh, "inlet", "outlet"; translate=(2, 0, 0))   # 1-based ids
```

```r
p <- mio_match_periodic_nodes(mesh, "inlet", "outlet", translate = c(2, 0, 0))
```

```js
const p = m.matchPeriodicNodes(mesh, 'inlet', 'outlet', [1,0,0,2, 0,1,0,0, 0,0,1,0, 0,0,0,1]);
```

Fortran, Julia and R ids are 1-based, C, Python and JavaScript 0-based. The MCP tool is `periodic`.
