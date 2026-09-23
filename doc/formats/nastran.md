# Nastran and OptiStruct (`.bdf`, `.fem`, `.nas`)

The [MSC/NX Nastran](https://help.autodesk.com/view/NSTRN/2019/ENU/?guid=GUID-42B54ACB-FBE3-47CA-B8FE-475E7AD91A00) bulk-data format — card entries (`GRID`, `CTRIA3`, `CTETRA`, `CHEXA`, …) in small-field, large-field or free (comma-separated) layout — including the [Altair OptiStruct](https://help.altair.com/hwsolvers/os/index.htm) dialect that HyperMesh writes, whose component names and sets become named regions.

| | |
|---|---|
| **Format name** | `nastran` |
| **Extensions** | `.bdf`, `.fem`, `.nas` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | — |

## Reading & writing

```python
import meshioplusplus

mesh = meshioplusplus.read("model.fem")
print([(r.name, r.kind, r.tag) for r in mesh.regions])  # HyperMesh components, SETs

meshioplusplus.nastran.write("out.bdf", mesh,
    point_format="fixed-large",   # "fixed-small", "fixed-large", or "free"
    cell_format="fixed-small",
)
```

- **`point_format`** / **`cell_format`** — the field layout for `GRID` and element cards. The C++ core writes the default `fixed-large`/`fixed-small` pair; any other choice uses the Python writer.

Both engines read any bulk-data deck and give the same mesh: points, cells, `nastran:ref` data and regions. Before v16.1.0 the C++ reader accepted only files the C++ writer had produced, so every real deck went through Python and the native CLI, C API and WASM could not read one.

## File structure

Everything before a line beginning `BEGIN BULK` — executive and case control, OptiStruct I/O options such as `OUTPUT,HM,ALL` — is skipped; `ENDDATA` ends the deck. Three field layouts are detected per line:

- **Small field (fixed)**: 10 fields of 8 characters; the tenth is the continuation marker.
- **Large field (fixed)**: 8 + 4×16 + 8 characters; the keyword carries a trailing `*` (`GRID*`, `CTETRA*`) and continuation lines start with `*`. The two 8-column halves of a 16-column field are re-merged before parsing.
- **Free field**: comma-separated, detected by a `,` on the line.

A card continues onto the next line when that line starts with `+` or `*`, or *implicitly* when the tenth field of one line and the first field of the next are both blank. Comment lines (`$`, `//`, `#`) may sit anywhere, including inside a continued card.

Real fields accept Nastran's compressed exponent (`1.5+1`, `.7E1`, `-2.45-16`); a blank real is 0.0. A malformed field raises a `ReadError` naming the card (`Nastran: invalid real field 'x' in a GRID card`).

`GRID`/`GRID*`: `[id, cp, x, y, z]`. Element cards: `[id, pid, nodes…]`. Each element card reads a fixed number of node fields and ignores whatever follows them — THETA/MCID and ZOFFS on shells, thicknesses on `CTRIA6`/`CQUAD8`, the orientation vector or grid of `CBAR`/`CBEAM`/`CBUSH`/`CBUSH1D`/`CGAP`. A missing node field is an error. The solids `CTETRA`, `CPYRA`/`CPYRAM`, `CPENTA` and `CHEXA` take their linear or quadratic form from the number of node fields given (4/10, 5/13, 6/15, 8/20); any other count is an error.

## Cell types & node ordering

| Nastran | meshio++ | Nastran | meshio++ |
|---|---|---|---|
| `CBEAM`, `CBUSH`, `CBUSH1D`, `CROD`, `CGAP`, `CBAR` | `line` | `CTETRA` (4 / 10 nodes) | `tetra` / `tetra10` |
| `CTRIAR`, `CTRIA3` | `triangle` | `CPYRAM`, `CPYRA` (5 / 13) | `pyramid` / `pyramid13` |
| `CTRAX6`, `CTRIAX6`, `CTRIA6` | `triangle6` | `CPENTA` (6 / 15) | `wedge` / `wedge15` |
| `CQUADR`, `CSHEAR`, `CQUAD4` | `quad` | `CHEXA` (8 / 20) | `hexahedron` / `hexahedron20` |
| `CQUAD8` | `quad8` | `CELAS1` (first grid) | `vertex` |
| `CQUAD9` | `quad9` | | |

The writer emits `CBAR`, `CTRIA3`, `CTRIA6`, `CQUAD4`, `CQUAD8`, `CQUAD9`, `CTETRA`, `CPYRA`, `CPENTA`, `CHEXA` and `CELAS1`, the solids with 10/13/15/20 grids for their quadratic forms. Before v16.1.0 both writers emitted the internal names `CTETRA_`, `CPYRA_`, `CPENTA_` and `CHEXA_` for those, which no solver accepts.

Node-order permutations, meshio++ slot `j` ← Nastran slot `P[j]`:

| type | permutation |
|---|---|
| `triangle6` from `CTRIAX6`/`CTRAX6` (corner, mid, corner, …) | `[0,2,4,1,3,5]`; `CTRIA6` already uses meshio++'s order |
| `hexahedron20` (vertical mid-edges before the top ring) | `[0,…,11,16,17,18,19,12,13,14,15]`, an involution |
| `wedge15` (likewise) | `[0,…,8,12,13,14,9,10,11]`, an involution |

## OptiStruct and HyperMesh

HyperMesh keeps what OptiStruct does not need in comment cards, which a plain Nastran reader discards. meshio++ reads:

- **Components.** `$HMMOVE <id>` is followed by `$` lines of element ids in 8-column fields, `a THRU b` for a range; `$HMNAME COMP <id>"name"` names the component and may come **after** the elements, so names are resolved at the end. Each component becomes a `cell` region named after it (or `component_<id>`), tagged with its id, with the dimension of its cells when they share one (else −1). A component that is named but has no members stays as an empty region.
- **Components by property.** HyperMesh writes `$HMMOVE` only for elements whose PID does not already say where they belong: `$HMNAME COMP <id>"name" <pid> "type" …` records the component's property after the name, and every element no `$HMMOVE` lists joins the component whose recorded property is its PID. Without a recorded property nothing is inferred.
- **Sets.** An OptiStruct `SET,<id>,GRID|ELEM,LIST,<ids>` card (ids and `THRU` ranges) becomes a `point` or `cell` region tagged with the set id, named by `$HMSET <id> <type> "name"` (or `set_<id>`). Other set types and subtypes are skipped with a warning.

Ids a component or set lists that name no grid or element are dropped with one warning; ids inside a `THRU` range pick only the ids that exist. `$HWCOLOR`, `$HMNAME PROP|MAT|LOADCOL|PLYS|LOADSTEP` and `$HMSETTYPE` are ignored.

**Skipped cards.** Cards that carry nothing meshio++ keeps are skipped without an error. Properties, materials, loads, constraints, coordinate systems, tables and solution parameters (`P…`, `MAT…`, `SPC…`, `MPC…`, `FORCE…`, `MOMENT…`, `LOAD…`, `TEMP…`, `GRAV`, `RFORCE`, `ACCEL…`, `CORD…`, `TABLE…`, `EIGR…`, `EIGC`, `NLPARM`, `TSTEP…`, `FREQ…`, `SUPORT`, `DAREA`, `DLOAD`, `RLOAD…`, `TLOAD…`, `SPOINT`, `ASET`, `OMIT`, `INCLUDE`) go quietly. Everything else — optimization (`DESVAR`, `DRESP1`, `DTPL`, …), contact (`CONTACT`, `TIE`, …), element cards meshio++ does not map (`CONM2`, `RBE2`, …) — is named, with counts, in a single warning at the end of the read:

```text
Nastran: skipped 6 card(s) meshio++ does not read: CONM2 (1), CONTACT (1), DESVAR (1), DRESP1 (2), TIE (1)
```

**Writing regions back.** Cell regions that are pairwise disjoint are written as HyperMesh components, before `ENDDATA`:

```text
$
$HMMOVE        5
$              1THRU          17
$HMNAME COMP                   5"a"
```

These are comments, so no solver sees them, and both engines write the same bytes. A region keeps its `tag` as the component id when every written region has a positive, unique one; otherwise they are numbered 1, 2, …. Point and side regions, and a cell region that overlaps one written before it, are dropped with a warning; regions are taken in `(kind, name, dim, tag)` order.

## Data mapping

- `point_data["nastran:ref"]` — the `GRID` card's CP field, when any card fills it; a blank field reads as 0 (the basic system).
- `cell_data["nastran:ref"]` — the element card's PID field, when any card fills it; a blank field reads as 0. The writer writes 0 as a blank field.
- `mesh.regions` — HyperMesh components and OptiStruct sets (above).
- `mesh.points_id` / `mesh.cells_id` — attributes (not data-dict entries) holding the original `GRID`/element ids; set by the Python reader only.

## Quirks & limitations

- The **16-character float encoding** of `GRID*` coordinates differs between the engines: the Python writer uses `np.format_float_scientific(precision=11)` with `E`; the C++ writer takes the shortest string (precision 0 to 11) that round-trips through `strtod`. Both fit the field; the text is not identical.
- The C++ writer's first line is the comment `$ meshioplusplus-cpp-nastran`, then the provenance tag; the Python writer writes only the tag. The sentinel is vestigial: the C++ reader of releases before 16.1 accepted only files carrying it.
- `CBAR`/`CBEAM`/`CBUSH`/`CBUSH1D`/`CGAP` lose their orientation vector or third grid on read, and the writer writes 2-noded `CBAR`s without one.
- Coordinates are read raw: a `GRID`'s CP system is kept in `nastran:ref`, not applied.
- A deck without `BEGIN BULK` (a bare include file) is refused.
- Points are promoted to 3D on write if given 2D, with a warning (Python writer).

## Notes

- `tests/python/meshes/nastran/cylinder.fem` and `cylinder_cells_first.fem` are HyperMesh 2017.3 exports matching `tests/python/meshes/med/cylinder.med`; their single component `misc1` holds every element. `optistruct_mixed.fem` is generated by `tools/gen_optistruct_fixture.py` (every layout, quadratic solids, components by `$HMMOVE` and by property, sets, skipped cards) and `composite_plate_2022.fem` is a real OptiStruct deck from [pyNastran](https://github.com/SteveDoyle2/pyNastran) (BSD-3-Clause). `tests/python/test_nastran.py` checks both engines against all four.
- The C++ reader uses `detail/keyword_card.hpp`'s `card_to_int`/`card_to_real` for its fields.
