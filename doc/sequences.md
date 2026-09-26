# Sequences: multi-file and transient datasets

Since v9.12.0 meshio++ can treat a **set of files** — or the steps inside one multi-step file — as one ordered logical dataset. That is how transient solver output actually arrives (`out_0000.vtu … out_0500.vtu`), and how most of the 41 formats have to express time, since only a minority carry several steps natively.

```bash
# fan-in: N single-step files -> one multi-step XDMF (quote the glob!)
meshioplusplus convert 'out_*.vtu' series.xdmf

# ... or one transient VTKHDF: the geometry once, one appended step per file
meshioplusplus convert 'out_*.vtu' series.vtkhdf

# fan-out: one multi-step file -> one file per step
meshioplusplus convert series.xdmf 'step_{step}.vtu'

# a whole transient post-processing run, chain applied per step
meshioplusplus pipeline transient.json
```

```python
import meshioplusplus as mp

for time, mesh in mp.read_sequence("out_*.vtu"):   # lazy: one mesh at a time
    print(time, mesh.points.shape)

mp.write_sequence("series.xdmf", mp.read_sequence("out_*.vtu"))
```

## It is a driver, not an operation

Everything here reads and writes through the existing format registry and runs operation chains through the existing **typed pipeline layer**. `run_pipeline_steps` remains the single owner of the step dispatch — there is deliberately no second `if (op == ...)` chain — so a sequence document and the browser viewer's `convertSurfaceOps` still cannot drift apart. Nothing here adds a mesh operation, a file format, or a dependency.

## The three shapes

| input expands to | output has `{step}`/`{index}` | mode |
|---|---|---|
| 1 file, 1 step | no | `sequence` (degenerate — identical to a plain single-file run) |
| 1 file, >1 step | yes | `fan-out` |
| 1 file, >1 step | no | **error** — refuses to truncate |
| >1 files | no | `fan-in` |
| >1 files | yes | `sequence` (N→N) |
| 1 file, 1 step | yes | `sequence` with one entry |

The mode is **inferred**. An explicit `Mode` (`"sequence"`, `"fan-in"`, `"fan-out"`) never *changes* the run: it **asserts** the inference and errors naming both on a mismatch. That is worth having because the inference depends on how many files a glob happened to match — a pattern matching exactly one file would otherwise quietly take the single-file path, and someone who wrote `"Mode": "fan-in"` wants that to fail.

![Fan-in, fan-out and N-to-N sequences, and the multi-step-to-single-file case that is refused rather than truncated](/diagrams/sequences_shapes.svg)

## Ordering is natural-numeric

A sequence has a defined order and it is **not lexicographic**: `out_10.vtu` must follow `out_9.vtu`, which a plain sort gets backwards. The rule is a documented contract, not an implementation detail:

1. Both strings split into maximal runs of digits and maximal runs of non-digits, compared run by run.
2. Two non-digit runs compare byte by byte as **unsigned** bytes — plain `char` signedness is implementation-defined, and would order the same UTF-8 paths differently on ARM than on x86.
3. Two digit runs compare **numerically**: leading zeros stripped, a shorter stripped run is less, equal lengths compare lexicographically. This is done on the digits themselves and never through `stoull`/`int()`, so a 40-digit hash-named file cannot overflow.
4. A digit run at the same position as a non-digit run sorts **first**.
5. If every run compares equal, the **unstripped** strings compare byte-wise, so `out_1` < `out_01` deterministically.
6. The **whole path** is compared, not the basename, so `a/out_9` < `a/out_10` < `b/out_1`.

Rule 5 is not cosmetic: without it `out_1` and `out_01` are mutually "not less" yet not equivalent, which makes the comparator not a strict weak ordering and `std::sort` undefined behaviour on any directory mixing padded and unpadded names. `tests/cpp/test_sequence.cpp` brute-forces all four axioms over a table containing exactly that case.

An **explicit path list** is a stated order and is *not* re-sorted; pass `sort=True` (Python) / `opts.sort` (C) to sort it anyway. A pattern is always sorted, since a directory listing has no meaningful order of its own.

## The pattern language

Exactly `*` (any run, possibly empty) and `?` (exactly one character). **No `**`, no `[set]`, no brace expansion**, and the directory component of a pattern is taken **literally** — a `*` there is an error naming the restriction rather than a recursive walk nobody asked for.

This is deliberately narrower than POSIX `glob(3)` *and* than Python's `glob`/`fnmatch`. Matching lives in the C++ core, so the native CLI, the C API, Fortran, Julia and R all get it; the Python twin re-implements the same restricted matcher rather than delegating to `fnmatch`, precisely so the two cannot accept different things. `[abc]` is three literal characters here.

### Shell quoting

**Quote the pattern.** A shell expands `out_*.vtu` before the CLI ever sees it, so `convert out_*.vtu out.xdmf` reaches the CLI as a dozen positional arguments and fails on the argument count. Two forms work:

```bash
meshioplusplus convert 'out_*.vtu' out.xdmf                       # quoted
meshioplusplus convert a.vtu --input b.vtu --input c.vtu out.xdmf  # pre-expanded
```

`--input` appends to the positional `infile` and is repeatable — it is what to reach for when your shell (or your build system) has already expanded the glob. Both forms produce the same dataset.

`--no-sequence` forces the single-file path, for the rare filename that genuinely contains `*` or `{step}` (legal on POSIX). `--sequence` forces the sequence path.

## Time values

Every entry carries a time and a record of **which source it came from**, because "the file said 0.25" and "nothing said anything, so this is position 3" are different facts and a user plotting against it needs to know which they have. The precedence:

| priority | source | reported as | where it comes from |
|---|---|---|---|
| 1 | explicit | `explicit` | the caller's own list (`Times` / `times=`) |
| 2 | file | `file` | a multi-step file's own step time, or a single-step file's `field_data["meshio:time"]` |
| 3 | filename | `filename` | the **last** maximal digit run of the stem, so `run17/out_0042.vtu` gives 42 |
| 4 | index | `index` | the integer position — the fallback, with a warning carrying the count |

`time_from` (`"auto"` by default) restricts this: `"file"`, `"filename"` and `"index"` each pin one source, falling back to the index when it says nothing.

An explicit list whose length does not match the entry count is an error naming both counts.

### `field_data["meshio:time"]`

A length-1 Float64 array holding the one time a mesh is a snapshot at. It generalizes the `exodus:time` convention the Exodus reader/writer already round-trips, and is colon-namespaced like `partition:part` because a bare `"time"` would collide with a solver's own field far too often.

**It only survives where the target format carries `field_data` at all**, which today means Exodus, Gmsh, MDPA and (since v15.0.0) VTU and VTP — VTK and most others carry none in either direction. A fan-out to `out_{step}.vtu` therefore records each step's time in the file (as `<FieldData>`, where ParaView also reads it as the time step), but a fan-in over the `{step}` glob still takes the number in the filename, which outranks a time found in the file, and the entry reports `filename` so you can see that. To close such a round trip on time, supply the times on the way back in:

```python
mp.write_sequence("series.xdmf", mp.read_sequence("step_*.vtu"), )   # times from filenames
mp.write_sequence("series.xdmf",
                  mp.read_sequence("step_*.vtu", times=[0.0, 0.1, 0.2]))
```

A fan-in to XDMF needs no help: the series records its own step times, and a later `sequence_entries` on it reports `file`.

## `{step}` and `{index}`

`{index}` is the plain decimal index. `{step}` is that index zero-padded to `max(4, digits(count - 1))`, so a 12-step run writes `out_0000 … out_0011` and a 20000-step run widens to `out_00005` rather than sorting wrongly in a directory listing.

Expansion is **substring replacement**, not `str.format`: an unrelated `{` in the path is a literal. This matches the native CLI's existing `{key}`/`{part}` helpers for `split` and `partition`, and differs from the *Python* CLI's `str.format`-based versions of those, which raise on a stray brace. That asymmetry is pre-existing; these tokens inherit the C++ side's semantics rather than introducing a third convention. (The three expanders are deliberately kept separate: `split`'s and `partition`'s output must stay byte-identical, and because the two languages differ a "shared" helper would be two helpers anyway.)

## Which formats carry time

Two different questions, answered two different ways:

**Reading** — how many steps does this file have? — is **registry-derived**: `read_metadata(...).time_values`. A format whose metadata reader fills no time values has one step, which is the truthful answer for every format that cannot express time. Today **XDMF**, **Exodus**, **GiD**, **MED**, **CGNS**, **Tecplot**, **Gmsh** (4.1 only), **EnSight**, **OpenFOAM**, **VTKHDF**, **PVD**, **CalculiX `.frd`**, **UNV**, **MSC Nastran HDF5**, **FEBio `.xplt`**, **Ansys `.rst`** (and its full-rotor reading `ansys_rst_cyclic`), **Marc `.t19`**, **LS-DYNA `d3plot`**, **Nastran OP2**, **OpenRadioss animation files** and **USD** report real counts. USD's is the one that is not cheap: it has no native metadata path, so the summary comes from a full read whose reader attaches the stage's authored sample times — correct, and stated here rather than implied. **MED**, **CGNS**, **Tecplot**, **Gmsh** and **EnSight** joined in v11.3.0 (roadmap §1 tier B1): `read_med_metadata` is a native path over `ENS_MAA`/`MAI` attributes and every `CHA/<field>/<step>`'s `PDT`, never a full read — the sorted, deduplicated union of every field's own step times, since a summary reports one timeline per file rather than one per field. `read_cgns_metadata` is the CGNS twin, over `Zone_t`'s own dimension triple and `BaseIterativeData_t/TimeValues`; the full CGNS reader now honours `time_step` too, picking the one `FlowSolution_t` the resolved step's `ZoneIterativeData_t/FlowSolutionPointers` names instead of reading (and silently letting a name collision between) every `FlowSolution_t` in the zone. `read_tecplot_metadata` scans every `ZONE` header's `SOLUTIONTIME`/`STRANDID` without decoding any zone's data body; the zones sharing the first zone's `STRANDID` (or, absent one, every zone with a `SOLUTIONTIME`) form the timeline `time_step` indexes into, sorted by `SOLUTIONTIME`, each distinct `SOLUTIONTIME` value giving one (possibly multi-zone) step. Several *non*-transient zones (no `SOLUTIONTIME` anywhere) are not a timeline, but are no longer truncated to the first either (v15.5.0, roadmap §1.1): the whole file becomes the single step, one cell block and one named [region](regions.md) per zone — see [Tecplot](formats/tecplot.md#multiple-zones). `read_gmsh_metadata` (4.1 only; 2.2 has no native metadata path and falls back regardless) reads the sorted union of every `$NodeData`/`$ElementData` section's own real-tag time value; the full reader resolves `time_step` against it and keeps exactly the sections matching that one time, replacing the old "first section per name wins" rule the default step (0) still reproduces unchanged. `read_ensight_metadata` reads a `.case` file's own `TIME` `time values:` (shape still needs a full read, like Exodus); unlike the other four, EnSight's own convention is **one file per step**, so `time_step` templates every `VARIABLE` entry's `*`-bearing filename rather than picking a section out of one file — see [EnSight](formats/ensight.md#reading-transient-variables). This closed the roadmap's last in-file-timeline item; reading a `VARIABLE` at all was new; the format was geometry-only before.

**VTKHDF** joined in v14.0.0: `read_vtkhdf_metadata`'s `time_values` come straight from `Steps/Values` — a native path that never touches the geometry, and for a polyhedron-free `UnstructuredGrid` the whole summary is native too (`fell_back_to_full_read = false`); `time_step` reads only that step's slice of each flat array through the `Steps` offset tables. See [VTKHDF](formats/vtkhdf.md#time).

**PVD** joined in v15.0.0: a `.pvd`'s steps are the distinct `timestep=` values of its index, so `read_pvd_metadata` reports `time_values` from the index alone — no piece is opened for them — and summarizes step 0's pieces for the rest; `time_step` selects the step and `piece` one `part` of it, and the chosen step's time is attached as `field_data["meshio:time"]`. A `.pvd` whose entries are `.pvtu` files is the ordinary layout of a partitioned transient run. See [PVD](formats/pvd.md).

**CalculiX `.frd`** joined in v15.3.0: its steps are the `100C` increments of the file, told apart by their frame id as well as their value (the modes of a frequency step can share a frequency). `read_frd_metadata` has no header-only path — it reads the file in full, like EnSight — so `fell_back_to_full_read` is true; `time_step` then parses only the chosen increment's result blocks, and the step's value, step counter and analysis type are attached as `field_data["meshio:time"]`, `frd:step` and `frd:analysis`. The format is read-only, so it is never a fan-in target. See [CalculiX results](formats/frd.md#results-and-steps).

**UNV** (`.unv`/`.uff`) joined in v15.6.0 (roadmap §1.1): its steps are the analysis steps of its results — every 2414/55/56 block is keyed by its analysis type, its step, mode or frequency number and its time, frequency or eigenvalue, blocks with the same key merge, and each abscissa value of the dataset-58 functions of a test-lab file (FRFs, coherences, spectra) is one step as well. `read_unv_metadata` reads the whole file, like `.frd`, and a step's value, analysis type and number are attached as `field_data["meshio:time"]`, `unv:analysis` and `unv:step`. The writer writes one step per file. See [UNV](formats/unv.md#results-and-steps).

**FEBio `.xplt`** joined in v16.2.0: its steps are the states of the plot file, one per saved time. `read_xplt_metadata` reads each state's header (inflating it when the file is compressed), and `time_step` then reads only that state's data. The state's time, index and status are attached as `field_data["meshio:time"]`, `xplt:step` and `xplt:status`, and `convert run.xplt 'out_{step}.vtu'` fans a run out to one `.vtu` per state. A remeshed run (v16.10.0) writes a new mesh before the states that use it, and each state is read on its own mesh, so a sequence's steps can differ in points and cells; `read_xplt_metadata` describes the first mesh. The format is read-only, so it is never a fan-in target. See [FEBio plot files](formats/xplt.md).

**Ansys `.rst`/`.rth`** joined in v16.3.0: its steps are the result sets of the file (load steps, substeps or modes), in the order of its data-set index. `read_ansys_rst_metadata` reads the headers and the model and takes every set's time, or frequency in a modal analysis, from the time table, so `convert file.rst 'mode_{step}.vtu'` writes one file per mode. `time_step` then reads only that set's results. `ansys_rst_cyclic` (v16.8.0) has the same steps, each the full rotor. See [`ansys_rst`](./formats/ansys_rst.md).

**Marc `.t19`** joined in v16.8.0: its steps are the increments of a Marc formatted post file, between its `****` and `----` markers. `read_marc_t19_metadata` reads each increment's time from its increment blocks (a frequency or buckling factor in a modal, harmonic or buckling analysis), so `convert job.t19 'inc_{step}.vtu'` writes one file per increment, and `time_step` reads only that increment's results. See [Marc](./formats/marc.md#the-post-file).

**Femap `.neu`** joined in v16.5.0: its steps are the `450` output sets, in file order, with the set value (a time, a frequency, a load factor, whatever the analysis stored) as the step time; each set's `451`/`1051` output vectors are its data. `read_femap_metadata` has no header-only path (it reads the file whole), and `femap_time_values` returns the values alone. A model with no output set is one step.

**Abaqus `.fil`** joined in v16.7.0 (roadmap §1.3): its steps are the increments of the results file (each `2000` … `2001` record pair), in file order, with the increment's total time as the step time; `abaqus:step`, `abaqus:increment`, `abaqus:step_time` and `abaqus:procedure` say which analysis step and increment it is. `read_abaqus_fil_metadata` has no header-only path: it reads the file once and takes every increment's time. See [Abaqus `.fil`](./formats/abaqus_fil.md#steps-and-data).

**MSC Nastran HDF5** (`.h5`) joined in v15.7.0 (roadmap §1.1): its steps are the result domains (`/NASTRAN/RESULT/DOMAINS` rows: a subcase, and within it a mode, time or frequency) that the INDEX tables of its nodal and element results reference, in `DOMAINS` order. `read_nastran_h5_metadata` reads the model and the first domain, like `.frd` (`fell_back_to_full_read` is true); `time_step` then reads only the chosen domain's rows, through hyperslabs. A step's `TIME_FREQ_EIGR` is `field_data["meshio:time"]` as written (an eigenvalue for a mode, not a frequency), with `nastran:subcase`, `nastran:mode` and `nastran:analysis` beside it. The format is read-only, so it is never a fan-in target. See [MSC Nastran HDF5](formats/nastran_h5.md#results-and-steps).

**LS-DYNA `d3plot`** joined in v16.9.0: its steps are the states of the family, `d3plot` then `d3plot01`, `d3plot02`... in numeric order, as many per file as fit before the file's end mark. `read_lsdyna_d3plot_metadata` reads the control block and the geometry and then only the time word of each state; `time_step` reads only the chosen state. The state's time and index are `field_data["meshio:time"]` and `lsdyna:state`. A glob over the family (`run/d3plot*`) keeps the base file only, so a d3plot folder is one sample, not one per file. The format is read-only, so it is never a fan-in target. See [LS-DYNA d3plot](formats/lsdyna_d3plot.md#states-and-data).

**Nastran OP2** joined in v16.9.0: its steps are the distinct (subcase, analysis code, mode, time or frequency) of the result tables it reads, in the order the file first names them. `read_nastran_op2_metadata` reads the model and the result headers; `time_step` decodes only the chosen step's records. The step's time, frequency or eigenvalue (0 for a static subcase) is `field_data["meshio:time"]`, with `nastran:subcase`, `nastran:analysis` and `nastran:mode`. See [Nastran OP2](formats/nastran_op2.md#results-and-steps).

**OpenRadioss animation files** (`radioss_anim`) joined in v16.11.0: a run writes one file per output time (`crashA001`, `crashA002`, …), each a single state whose time is in the file. `read_radioss_anim_metadata` reports that one time (it reads the file whole), so a glob over the run, `convert 'crashA*' crash.vtkhdf`, takes every file's own time (`time_source` `file`) rather than the number in its name. The state's time is `field_data["meshio:time"]`. The format is read-only, so it is never a fan-in target. See [OpenRadioss animation files](formats/radioss_anim.md).

**LS-DYNA binout** (`lsdyna_binout`) joined in v16.12.0: its steps are the `nodout` outputs (else the first database's), each a point cloud of the `nodout` nodes; every other database's latest output at or before the step's time is field data. `read_lsdyna_binout_metadata` reads the whole archive (`fell_back_to_full_read` is true). See [LS-DYNA binout](formats/lsdyna_binout.md#steps-and-data).

**OpenRadioss time-history files** (`radioss_th`) joined in v16.12.0: its steps are the outputs of a `<run>T01` file, each field data only (no points): the time, the global, part and subset variables and one array per TH group. A run stopped mid-write leaves an incomplete last output, which is not a step. See [OpenRadioss time-history files](formats/radioss_th.md#steps-and-data).

**DOLFINx VTX** (`vtx`) joined in v16.13.0: its steps are the ADIOS2 steps of a `.bp` directory, each time its `step` value; a step without mesh variables (a `VTXMeshPolicy.reuse` file) takes the last mesh written. `read_vtx_metadata` reads the step times and the block sizes without reading points (`fell_back_to_full_read` is false). See [DOLFINx VTX](formats/vtx.md#mapping).

**Tecplot SZL** (`szplt`) joined in v16.13.0: its steps are `tecplot`'s, the distinct solution times, read through TecIO. See [Tecplot SZL](formats/szplt.md).

**OpenFOAM** joined in v11.4.0 (roadmap §1 tier B2): `read_openfoam_metadata`'s `time_values` is a real, cheap native path — a directory listing of the case root's numeric-named subdirectories that hold at least one field file, sorted ascending — but the rest of that summary still needs a full read (`fell_back_to_full_read = true`), like Exodus/EnSight, since a cheap point/cell count would otherwise re-derive the whole cell-reconstruction pipeline. `time_step` selects one time directory (an index into the sorted list, not the directory's own name); `arrays` selects which of its field files to attach. Unlike the five above, OpenFOAM's "steps" are directories of independent files rather than sections of one — closer to EnSight's own convention, except every step lives in a fully separate `<time>/` tree rather than a numbered filename. See [OpenFOAM](formats/openfoam.md#time-directory-fields).

**Writing** — can this format hold N steps? — has no file to probe, so it is a small owned predicate: **XDMF, GiD, USD, VTKHDF, PVD and Femap** (Femap since v16.17.0: a neutral file holds one mesh and any number of output sets, so a fan-in is the mesh once plus one output set per entry, and every step must have the first step's cells — see [Femap](formats/femap.md#time-series). PVD since v15.0.0: a collection is an index over one `.vtu` per step, so a fan-in writes one piece and one index line per entry and rewrites the index after every step, so a killed run still leaves a collection ParaView opens; `Encoding` (ASCII or binary pieces) is honoured, `Codec` and `FloatFormat` are rejected — see [PVD](formats/pvd.md). VTKHDF since v14.0.0: its `Steps` group is an offset table into flat arrays, so a fan-in is the geometry written once plus one appended step per entry — see [VTKHDF time series](vtkhdf_time_series.md)). A USD stage is the natural fan-in target of the three: points are time-sampled per step and topology re-authored only when it changes, so a remeshed series stays valid while a fixed-topology one stays compact. The anti-drift mechanism is a test rather than the table: a gtest iterates every registered writer and asserts the predicate agrees with whether a real two-step fan-in to that format actually succeeds, so a format that grows a series writer without updating the predicate turns CI red naming itself.

**Fan-out is the answer for everything else**, which is most of the format table.

> **Directory stores in a glob.** `pmsh` and `zarr` write a *directory* per step, and a DOLFINx VTX file is a `.bp` directory, which an ordinary file glob would skip. `read_sequence("out_*.pmsh")` and `read_sequence("run_*.bp")` match them anyway: a plan entry is a file, or a directory whose suffix is one of those three. Deliberately narrower than "any directory", so an unrelated subdirectory that happens to match a pattern is still ignored. The C++ `sequence_expand` keeps `.bp` directories (v16.13.0) and not the other two, which costs nothing there — `pmsh` and `zarr` are Python-only and unreadable from that surface. An Elmer mesh directory has no suffix, so both globs keep a matching directory by its files instead (`mesh.header`, or a `partitioning.N` child), since v16.17.0: `read_sequence("step_*")` over `step_1/`, `step_2/`… reads one Elmer mesh per step. An OpenFOAM case is not kept, since it holds several steps itself.

## Never a silent truncation

A multi-step input aimed at a single-step output is an **error naming the remedy**, not a quiet write of step 0:

```
meshio++: sequence: the input has 12 time steps but the output path 'one.vtu'
names a single file; add '{step}' to write one file per step, or select a
single step (--time-step on the CLI, Input.Options.TimeStep in a settings
document)
```

and a fan-in to a format that cannot hold a series:

```
meshio++: sequence: format 'vtu' cannot hold a multi-step series (only 'xdmf',
'gid', 'vtkhdf', 'pvd' and 'femap' can); write one file per step with an Output path containing
'{step}' instead
```

An explicit `--time-step` / `Input.Options.TimeStep` **is** a deliberate single-step selection and opts out of the check. Before v9.12.0, `convert series.xdmf out.vtu` silently wrote the first step.

## The streaming guarantee

> At most **one** `Mesh` is alive at any point inside a fan-in, a fan-out or a per-step run, and `read_sequence` yields at most one live mesh per iteration.

This is a **contract, not an optimization** — the whole feature exists so a 500-step dataset is traversable on a laptop. No implementation may buffer the sequence, including "just for sorting" or "just to compute the time range": expansion returns the *plan* (paths, step indices and times), never meshes.

It is pinned rather than asserted in prose. A gtest measures the peak through the `BufferAllocator` hook and requires it to be **O(1) in the step count** (the real assertion is that the 40-file peak matches the 20-file one, not that either is small); the Python suite uses weak references, so a regression names the retainer rather than merely being slow.

The C ABI expresses the same rule: `mio_sequence_read` hands back an **owned** mesh, deliberately unlike `mio_split_result_mesh`'s borrow, because a borrow would force the handle to cache every mesh it produced.

**What does not stream**: `list(read_sequence(...))` is your choice and your memory, and an operation's own internals are unchanged — a single step that does not fit in memory still does not fit.

## Holding a series as one value: `TimeSeries`

`read_sequence` is a **generator**: exhausted after one pass, unable to answer `len(...)`, unable to be indexed. That is deliberate — it is the surface for a single streaming traversal. But the C/Fortran/Julia/R sequence handles (`mio_sequence`, `type(mio_sequence)`, `Sequence`, `mio_sequence()`) already give **random access** by construction: they hold the entry *plan* and read step *i* on demand, any number of times, in any order. Python had no equivalent — until `TimeSeries`:

```python
series = mp.TimeSeries("out_*.vtu")
len(series)                  # 12, from the plan alone -- no reads
series.times                 # [0.0, 1.0, ..., 11.0], likewise
t0, mesh0 = series[0]        # exactly one read
t_last, mesh_last = series[-1]
for t, mesh in series:       # a fresh, independent pass every time
    ...
```

It closes the one item [`doc/roadmap.md`](roadmap.md) had left open for this feature: a caller that genuinely needs to *hold* a series — random access across steps, more than one pass — no longer has to reach for `list(...)`. It still honours the streaming invariant: only the plan is held, never a mesh, so `series[i]` performs exactly one read and nothing is cached between accesses — holding a 500-entry `TimeSeries` costs no more memory than holding its plan.

## Sequences in a settings document

The v9.11.0 pipeline schema, plus seven keys. See [the settings pipeline](pipeline.md) for everything else.

```jsonc
{
  "Version": 1,
  "Mode": "sequence",                 // optional: asserts the inferred shape
  "Input": {
    "Pattern": "run/out_*.vtu",       // or "Paths": [...], or "Path": "in.xdmf"
    "Times": [0.0, 0.1, 0.2],         // optional explicit times
    "TimeFrom": "auto",               // auto | file | filename | index
    "Format": "vtu",
    "Options": { "PointsOnly": false }
  },
  "Operations": [                     // applied to EVERY step
    { "Op": "Quality" },
    { "Op": "Gradient", "Array": "temperature" }
  ],
  "Output": { "Path": "post/out_{step}.vtu" },
  "Parallel": true,                   // Python driver only
  "Workers": 8
}
```

`Path`, `Pattern` and `Paths` are mutually exclusive; naming more than one is an error. Parsing stays strict everywhere: an unknown key, op or enum value errors naming it.

**Backward compatibility is absolute.** A document using none of the new keys and naming a plain output takes a physically unchanged path — the C++ engine literally delegates to `run_pipeline`. Conversely, `parse_pipeline_json` (the typed single-file parser) **rejects** a sequence key by name rather than ignoring it, because ignoring one would run a transient document as its first step.

### `Parallel` / `Workers`

Files are embarrassingly parallel at the driver level, so the **Python** driver can run them in a `ProcessPoolExecutor` (`Workers` defaults to one per core). Output is identical to the serial run, report included.

Two deliberate restrictions:

- **`Parallel` with a fan-in is an error**, not a silent serialization. The series writer is one stateful handle whose steps must be appended in order, and the obvious workaround — parallel reads feeding an ordered queue into one writer — would hold `Workers + 1` meshes and break the streaming guarantee. If that trade is ever worth making it should be an explicit, named opt-in with its memory cost stated.
- **The C++ engine runs serially** and records a warning when a document asks for parallelism. Every operation already parallelizes internally, so a step-level parallel region nested over them would oversubscribe on every backend.

## Per-surface entry points

| Surface | Call |
|---|---|
| Python CLI | `meshioplusplus convert 'in_*.vtu' out.xdmf`, `… in.xdmf 'out_{step}.vtu'`, `meshioplusplus pipeline settings.json` |
| Native CLI | the same words |
| Python | `read_sequence`, `write_sequence`, `sequence_entries`, `run_sequence_pipeline`, `TimeSeries` (and `run_pipeline`, which routes here) |
| C | `mio_sequence_open`/`_open_list`/`_count`/`_path`/`_step`/`_time`/`_time_source`/`_read`/`_free`, `mio_sequence_to_timeseries`/`_ex`, `mio_timeseries_to_sequence`, `mio_sequence_pipeline_run_file`/`_json` |
| Fortran | `type(mio_sequence)` with `%open`/`%count`/`%path`/`%time`/`%read_step`/`%to_timeseries`/`%free`, plus `mio_timeseries_to_sequence` |
| Julia | `Sequence`, `read_step`, `to_timeseries`, `timeseries_to_sequence`, `run_sequence_file`/`_json` |
| R | `mio_sequence()`, `mio_sequence_read()`, `mio_sequence_to_timeseries()`, `mio_timeseries_to_sequence()`, … |
| MCP | the `sequence` tool |
| WASM | `sequenceEntries`, `openSequence` (stateful, per-step reads), `sequenceToTimeseries`, `timeseriesToSequence`, and `runPipeline` (which routes) — over MEMFS paths |

### WASM

The browser surface has these too, over MEMFS paths — the same filesystem `convert` and `runPipeline` already work on:

```js
m.sequenceEntries('/seq/out_*.vtu');                      // the ordered plan
m.sequenceToTimeseries('/seq/out_*.vtu', '/seq/s.xdmf');  // fan-in
m.timeseriesToSequence('/seq/s.xdmf', '/seq/b_{step}.vtu');  // fan-out -> paths

// Stateful, lazy per-step reads (mirrors the C API's mio_sequence_* handle):
const seq = m.openSequence('/seq/out_*.vtu');
seq.read(3, { pointsOnly: true });   // readMeshSelective's options, per step
seq.close();
```

`Parallel` is accepted and **ignored with a warning** there: it is a Python-driver process pool, and a wasm module has no processes to pool. See [the WASM docs](wasm.md#sequences-transient--multi-file-datasets).

## Worked examples

### Convert a directory of steps into one time series

```bash
meshioplusplus convert 'results/out_*.vtu' results/series.xdmf
```

Reads the twelve `out_0.vtu … out_11.vtu` in natural-numeric order (so `out_10` follows `out_9`), takes each time from the trailing digit run of its filename, and streams them into one temporal XDMF collection — one mesh in memory at a time, whatever the step count.

### A transient post-processing run

`transient.json`:

```json
{
  "Version": 1,
  "Input": { "Pattern": "raw/out_*.vtu" },
  "Operations": [
    { "Op": "Quality" },
    { "Op": "Gradient", "Array": "temperature", "Operator": "gradient" },
    { "Op": "Clean" }
  ],
  "Output": { "Path": "post/out_{step}.vtu" },
  "Parallel": true
}
```

```bash
meshioplusplus pipeline transient.json
```

Every step is read, run through the same three-operation chain, and written to its own file — the composition that turns the pipeline into a batch post-processor. The report carries one entry per (step, op).

### Explode a multi-step file for a tool that cannot read one

```bash
meshioplusplus convert simulation.xdmf 'steps/frame_{step}.vtu'
```

Writes `steps/frame_0000.vtu …`, each carrying its own `field_data["meshio:time"]` where the target format can hold it.

## See also

- [The settings pipeline](pipeline.md) — the operation chain a sequence applies per step.
- [XDMF time series](xdmf_time_series.md) — the stateful writer a fan-in drives.
- [VTKHDF time series](vtkhdf_time_series.md) — the other one: `Steps` is an offset table, so a fan-in is geometry × 1, fields × N.
- [Selective reads](selective_read.md) — `time_step`, and what `read_metadata` can answer without a full read.
- [CLI reference](cli.md#meshioplusplus-convert).
- [Dataset manifests](datasets.md) — cataloguing many sequences (each entry *is* a `TimeSeries` plan) with splits/tags for ML training.
