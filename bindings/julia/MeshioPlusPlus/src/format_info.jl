# Format side channels: what a format carries that a mesh cannot hold (MDPA's
# tables, geometries, Mesh blocks and constraints), kept by `read_with_info` in
# an opaque handle for `write_with_info` to put back. The C API's
# `mio_format_info` (v16.27.0).

const MIO_MDPA_PROPERTIES = Int32(0)
const MIO_MDPA_ENTITY_NAMES = Int32(1)
const MIO_MDPA_SKIPPED = Int32(2)
const MIO_MDPA_MODEL_PART_DATA = Int32(3)
const MIO_MDPA_TABLES = Int32(4)
const MIO_MDPA_GEOMETRIES = Int32(5)
const MIO_MDPA_MESH_BLOCKS = Int32(6)
const MIO_MDPA_SUBMODELPARTS = Int32(7)
const MIO_MDPA_RAW_BLOCKS = Int32(8)
const MIO_GMSH_BOUNDING_ENTITIES = Int32(0)
const MIO_GMSH_PERIODIC = Int32(1)
const MIO_NURBS_PATCHES = Int32(0)
const MIO_NURBS_FIELDS = Int32(1)
const MIO_NURBS_TOPOLOGY = Int32(2)
const MIO_NURBS_CELL_MAP = Int32(3)
const MIO_BREP_FILE = Int32(0)
const MIO_BREP_BREPS = Int32(1)
const MIO_BREP_PATCHES = Int32(2)
const MIO_BREP_FACES = Int32(3)
const MIO_BREP_LOOPS = Int32(4)
const MIO_BREP_TRIMS = Int32(5)
const MIO_BREP_EDGES = Int32(6)
const MIO_BREP_VERTICES = Int32(7)
const MIO_BREP_REFS = Int32(8)

"""
    FormatInfo

A format's side channel, returned by [`read_with_info`](@ref). Released by a
finalizer; [`close`](@ref) does it deterministically and is idempotent. Inspect
an MDPA one with [`mdpa_info`](@ref), a Gmsh one with [`gmsh_info`](@ref), an
MFEM NURBS one with [`nurbs_info`](@ref), or a B-rep one (Kratos CAD JSON, STEP,
IGES) with [`brep_info`](@ref).
"""
mutable struct FormatInfo
    ptr::Ptr{Cvoid}

    function FormatInfo(ptr::Ptr{Cvoid})
        info = new(ptr)
        finalizer(info) do x
            if x.ptr != C_NULL
                ccall(_sym(:mio_format_info_free), Cvoid, (Ptr{Cvoid},), x.ptr)
                x.ptr = C_NULL
            end
        end
        info
    end
end

function Base.close(info::FormatInfo)
    if info.ptr != C_NULL
        ccall(_sym(:mio_format_info_free), Cvoid, (Ptr{Cvoid},), info.ptr)
        info.ptr = C_NULL
    end
    nothing
end

function _handle(info::FormatInfo)
    info.ptr == C_NULL && throw(MeshioError(MIO_ERR_INVALID_ARG, "FormatInfo is closed"))
    info.ptr
end

"""
    read_with_info(path; format="", options=nothing) -> (Mesh, Union{FormatInfo,Nothing})

Read a mesh, keeping what it cannot hold. For a format with a side channel
(currently `"mdpa"`, `"gmsh"` and `"mfem"`) the second value is a [`FormatInfo`](@ref); for every other
format it is `nothing` and the read is exactly [`read`](@ref)'s.
"""
function read_with_info(path::AbstractString; format::AbstractString="",
                        options::Union{Nothing,ReadOptions}=nothing)
    info = Ref{Ptr{Cvoid}}(C_NULL)
    ptr = _with_read_opts(options === nothing ? ReadOptions() : options) do ref
        ccall(_sym(:mio_read_with_info), Ptr{Cvoid},
              (Cstring, Cstring, Ptr{_CReadOpts}, Ptr{Ptr{Cvoid}}), path, format, ref, info)
    end
    mesh = Mesh(_check_ptr(ptr))
    mesh, info[] == C_NULL ? nothing : FormatInfo(info[])
end

"""
    write_with_info(mesh, info, path; format="")

Write a mesh, restoring the side channel `info` came back with. `info ===
nothing` writes exactly as [`write`](@ref); a side channel of another format
throws.
"""
function write_with_info(m::Mesh, info::Union{FormatInfo,Nothing}, path::AbstractString;
                         format::AbstractString="")
    h = info === nothing ? C_NULL : _handle(info)
    _check(ccall(_sym(:mio_write_with_info), Cint, (Cstring, Ptr{Cvoid}, Cstring, Ptr{Cvoid}),
                 path, _handle(m), format, h))
end

"""
    format_name(info) -> String

The format the side channel belongs to (`"mdpa"`, `"gmsh"` or `"mfem"`).
"""
format_name(info::FormatInfo) =
    _getstring((buf, n) -> ccall(_sym(:mio_format_info_format), Int64,
                                 (Ptr{Cvoid}, Ptr{UInt8}, Int64), _handle(info), buf, n))

_mdpa_count(h, section) = Int(_check_count(ccall(_sym(:mio_mdpa_info_count), Int64,
                                                 (Ptr{Cvoid}, Int32), h, section), "mdpa count"))

_mdpa_string(h, section, i, field) =
    _getstring((buf, n) -> ccall(_sym(:mio_mdpa_info_string), Int64,
                                 (Ptr{Cvoid}, Int32, Int64, Int32, Ptr{UInt8}, Int64),
                                 h, section, Int64(i), Int32(field), buf, n))

function _mdpa_int(h, section, i, field)
    v = Ref{Int64}(0)
    _check(ccall(_sym(:mio_mdpa_info_int), Cint, (Ptr{Cvoid}, Int32, Int64, Int32, Ptr{Int64}),
                 h, section, Int64(i), Int32(field), v))
    v[]
end

"""Copy one borrowed array out of the handle, in C (row-major) order, reshaped
column-major: an `(rows, cols)` C array comes back as a `(cols, rows)` matrix."""
function _mdpa_copy(f)
    data = Ref{Ptr{Cvoid}}(C_NULL); dt = Ref{Cint}(0); nd = Ref{Int32}(0)
    shape = zeros(Int64, MIO_MAX_NDIM)
    _check(f(data, dt, nd, shape))
    T = _jltype(dt[])
    dims = Tuple(Int.(reverse(shape[1:nd[]])))
    total = prod(dims; init=1)
    total == 0 && return zeros(T, dims)
    copy(unsafe_wrap(Array, Ptr{T}(data[]), dims))
end

_mdpa_array(h, section, i, field) =
    _mdpa_copy((d, t, n, s) -> ccall(_sym(:mio_mdpa_info_array), Cint,
                                     (Ptr{Cvoid}, Int32, Int64, Int32, Ptr{Ptr{Cvoid}}, Ptr{Cint},
                                      Ptr{Int32}, Ptr{Int64}),
                                     h, section, Int64(i), Int32(field), d, t, n, s))

"""
    gmsh_info(info) -> NamedTuple

Copy of Gmsh `bounding_entities` (signed entity tags per cell block) and
`periodic` links `(entity, affine, node_pairs)`. Entity is `(dimension, slave tag,
master tag)` in file ids. Node pairs are 1-based point rows in a `(2, N)` matrix;
ordering/duplicates are kept. Mesh operations do not remap the side channel.
"""
function gmsh_info(info::FormatInfo)
    h = _handle(info)
    format_name(info) == "gmsh" ||
        throw(MeshioError(MIO_ERR_INVALID_ARG, "not a gmsh side channel"))
    n(section) = Int(_check_count(ccall(_sym(:mio_gmsh_info_count), Int64,
                                       (Ptr{Cvoid}, Int32), h, section), "gmsh count"))
    a(section, i, field) = _mdpa_copy((d, t, nd, s) -> ccall(_sym(:mio_gmsh_info_array), Cint,
        (Ptr{Cvoid}, Int32, Int64, Int32, Ptr{Ptr{Cvoid}}, Ptr{Cint}, Ptr{Int32}, Ptr{Int64}),
        h, section, Int64(i), Int32(field), d, t, nd, s))
    (bounding_entities=[a(MIO_GMSH_BOUNDING_ENTITIES, i, 0)
                        for i in 0:n(MIO_GMSH_BOUNDING_ENTITIES)-1],
     periodic=[(entity=Tuple(a(MIO_GMSH_PERIODIC, i, 0)),
                affine=a(MIO_GMSH_PERIODIC, i, 1),
                node_pairs=a(MIO_GMSH_PERIODIC, i, 2) .+ 1)
               for i in 0:n(MIO_GMSH_PERIODIC)-1])
end

"""
    nurbs_info(info) -> NamedTuple

Copy of an MFEM NURBS side channel (v16.40.0): the patches an `MFEM NURBS mesh`
was sampled from. `patches` holds `(degrees, knots, control_points, weights,
fields)` per patch, `control_points` a `(space_dimension, points)` matrix with
direction 1 fastest and each field a `(components, points)` matrix of control
values; then `field_names`, `topology` (`nothing` when the source states none;
else `(dimension, num_vertices, patch_vertices, patch_attributes,
boundary_vertices, boundary_attributes, edges)`, `edges` as `(knot vector,
vertex, vertex)` columns) and the cell map `cell_patch`, `cell_span` (a
`(3, cells)` matrix) and `boundary_cell_patch`. Vertices, patches and knot
spans are 1-based; attributes and knot-vector ids stay as the file spelled
them. [`write_with_info`](@ref) writes the patches back from the handle.
"""
function nurbs_info(info::FormatInfo)
    h = _handle(info)
    format_name(info) == "mfem" ||
        throw(MeshioError(MIO_ERR_INVALID_ARG, "not an mfem (NURBS) side channel"))
    n(section) = Int(_check_count(ccall(_sym(:mio_nurbs_info_count), Int64,
                                       (Ptr{Cvoid}, Int32), h, section), "nurbs count"))
    a(section, i, field) = _mdpa_copy((d, t, nd, s) -> ccall(_sym(:mio_nurbs_info_array), Cint,
        (Ptr{Cvoid}, Int32, Int64, Int32, Ptr{Ptr{Cvoid}}, Ptr{Cint}, Ptr{Int32}, Ptr{Int64}),
        h, section, Int64(i), Int32(field), d, t, nd, s))
    function ival(section, i, field)
        v = Ref{Int64}(0)
        _check(ccall(_sym(:mio_nurbs_info_int), Cint, (Ptr{Cvoid}, Int32, Int64, Int32, Ptr{Int64}),
                     h, section, Int64(i), Int32(field), v))
        v[]
    end
    names = [_getstring((buf, len) -> ccall(_sym(:mio_nurbs_info_string), Int64,
                                            (Ptr{Cvoid}, Int32, Int64, Int32, Ptr{UInt8}, Int64),
                                            h, MIO_NURBS_FIELDS, Int64(k), Int32(0), buf, len))
             for k in 0:n(MIO_NURBS_FIELDS)-1]
    patches = [(degrees=a(MIO_NURBS_PATCHES, p, 0),
                knots=[a(MIO_NURBS_PATCHES, p, 3 + d) for d in 0:ival(MIO_NURBS_PATCHES, p, 0)-1],
                control_points=a(MIO_NURBS_PATCHES, p, 1),
                weights=a(MIO_NURBS_PATCHES, p, 2),
                fields=[a(MIO_NURBS_PATCHES, p, 6 + k) for k in 0:length(names)-1])
               for p in 0:n(MIO_NURBS_PATCHES)-1]
    topology = if n(MIO_NURBS_TOPOLOGY) == 0
        nothing
    else
        edges = a(MIO_NURBS_TOPOLOGY, 0, 4)
        edges[2:3, :] .+= 1
        (dimension=ival(MIO_NURBS_TOPOLOGY, 0, 0), num_vertices=ival(MIO_NURBS_TOPOLOGY, 0, 1),
         patch_vertices=a(MIO_NURBS_TOPOLOGY, 0, 0) .+ 1,
         patch_attributes=a(MIO_NURBS_TOPOLOGY, 0, 1),
         boundary_vertices=a(MIO_NURBS_TOPOLOGY, 0, 2) .+ 1,
         boundary_attributes=a(MIO_NURBS_TOPOLOGY, 0, 3), edges=edges)
    end
    cells(field) = n(MIO_NURBS_CELL_MAP) == 0 ? Int64[] : a(MIO_NURBS_CELL_MAP, 0, field) .+ 1
    (patches=patches, field_names=names, topology=topology,
     cell_patch=cells(0), cell_span=cells(1), boundary_cell_patch=cells(2))
end

"""
    brep_info(info) -> NamedTuple

Copy of a B-rep side channel: a Kratos CAD JSON handle (`"cad_json"`), or a
STEP/IGES one read through OpenCASCADE. `patches` holds `(degrees, knots,
control_points, weights)` per face surface, `control_points` a `(3, points)`
matrix with direction 1 fastest and the knots in the full form; `faces` holds
`(id, name, patch, swapped_normal, trimmed, rational, reduced_knots,
control_point_ids, loops, embedded_edges, surface_extra, extra)`, each loop
`(type, trims, extra)` and each trim `(trim_index, direction, curve, extra)`;
`edges` holds `(id, name, curve, topology, extra)` (`curve` is `nothing` for an
edge without a 3-D curve) and `vertices` `(id, name, topology, extra)`, each
topology entry `(id, name, trim_index, relative_direction, local_coordinates,
extra)`; `breps` holds `(id, name, faces, edges, vertices, extra)`; and
`extra` is the file's other members. A curve is `(degree, knots, reduced_knots,
rational, control_points, weights, ids, active_range, extra)`. Indices (a
face's `patch`, a brep's lists) are 1-based; ids (-1 where the file names the
entity instead), trim indices and the `extra` JSON texts stay as the file
spelled them. [`write_with_info`](@ref) writes the B-rep back as CAD JSON.
"""
function brep_info(info::FormatInfo)
    h = _handle(info)
    format_name(info) in ("cad_json", "step", "iges") ||
        throw(MeshioError(MIO_ERR_INVALID_ARG, "not a B-rep side channel"))
    n(section) = Int(_check_count(ccall(_sym(:mio_brep_info_count), Int64,
                                       (Ptr{Cvoid}, Int32), h, section), "brep count"))
    a(section, i, field) = _mdpa_copy((d, t, nd, s) -> ccall(_sym(:mio_brep_info_array), Cint,
        (Ptr{Cvoid}, Int32, Int64, Int32, Ptr{Ptr{Cvoid}}, Ptr{Cint}, Ptr{Int32}, Ptr{Int64}),
        h, section, Int64(i), Int32(field), d, t, nd, s))
    function ival(section, i, field)
        v = Ref{Int64}(0)
        _check(ccall(_sym(:mio_brep_info_int), Cint, (Ptr{Cvoid}, Int32, Int64, Int32, Ptr{Int64}),
                     h, section, Int64(i), Int32(field), v))
        v[]
    end
    str(section, i, field) = _getstring((buf, len) -> ccall(_sym(:mio_brep_info_string), Int64,
                                                            (Ptr{Cvoid}, Int32, Int64, Int32,
                                                             Ptr{UInt8}, Int64),
                                                            h, section, Int64(i), Int32(field),
                                                            buf, len))
    curve(section, i) = (degree=Int(ival(section, i, 10)), knots=a(section, i, 12),
                         reduced_knots=ival(section, i, 12) != 0,
                         rational=ival(section, i, 11) != 0, control_points=a(section, i, 10),
                         weights=a(section, i, 11), ids=a(section, i, 13),
                         active_range=a(section, i, 14), extra=str(section, i, 10))
    trim(t) = (trim_index=ival(MIO_BREP_TRIMS, t, 2), direction=ival(MIO_BREP_TRIMS, t, 3) != 0,
               curve=curve(MIO_BREP_TRIMS, t), extra=str(MIO_BREP_TRIMS, t, 0))
    ref(r) = (id=ival(MIO_BREP_REFS, r, 0), name=str(MIO_BREP_REFS, r, 0),
              trim_index=ival(MIO_BREP_REFS, r, 1), relative_direction=ival(MIO_BREP_REFS, r, 2),
              local_coordinates=a(MIO_BREP_REFS, r, 0), extra=str(MIO_BREP_REFS, r, 1))
    span(first, count) = first:first+count-1
    patches = [(degrees=Int32[ival(MIO_BREP_PATCHES, p, 0), ival(MIO_BREP_PATCHES, p, 1)],
                knots=[a(MIO_BREP_PATCHES, p, 2), a(MIO_BREP_PATCHES, p, 3)],
                control_points=a(MIO_BREP_PATCHES, p, 0), weights=a(MIO_BREP_PATCHES, p, 1))
               for p in 0:n(MIO_BREP_PATCHES)-1]
    faces = map(0:n(MIO_BREP_FACES)-1) do f
        iv(k) = ival(MIO_BREP_FACES, f, k)
        loops = [(type=str(MIO_BREP_LOOPS, l, 0),
                  trims=[trim(t) for t in span(ival(MIO_BREP_LOOPS, l, 1), ival(MIO_BREP_LOOPS, l, 2))],
                  extra=str(MIO_BREP_LOOPS, l, 1))
                 for l in span(iv(6), iv(7))]
        embedded = [trim(t) for t in span(iv(8), iv(9)) if ival(MIO_BREP_TRIMS, t, 1) < 0]
        (id=iv(0), name=str(MIO_BREP_FACES, f, 0), patch=Int(iv(1)) + 1,
         swapped_normal=iv(2) != 0, trimmed=iv(3) != 0, rational=iv(4) != 0,
         reduced_knots=iv(5) != 0, control_point_ids=a(MIO_BREP_FACES, f, 0), loops=loops,
         embedded_edges=embedded, surface_extra=str(MIO_BREP_FACES, f, 2),
         extra=str(MIO_BREP_FACES, f, 1))
    end
    edges = [(id=ival(MIO_BREP_EDGES, e, 0), name=str(MIO_BREP_EDGES, e, 0),
              curve=ival(MIO_BREP_EDGES, e, 1) != 0 ? curve(MIO_BREP_EDGES, e) : nothing,
              topology=[ref(r) for r in span(ival(MIO_BREP_EDGES, e, 2), ival(MIO_BREP_EDGES, e, 3))],
              extra=str(MIO_BREP_EDGES, e, 1))
             for e in 0:n(MIO_BREP_EDGES)-1]
    vertices = [(id=ival(MIO_BREP_VERTICES, v, 0), name=str(MIO_BREP_VERTICES, v, 0),
                 topology=[ref(r) for r in span(ival(MIO_BREP_VERTICES, v, 1),
                                                ival(MIO_BREP_VERTICES, v, 2))],
                 extra=str(MIO_BREP_VERTICES, v, 1))
                for v in 0:n(MIO_BREP_VERTICES)-1]
    breps = [(id=ival(MIO_BREP_BREPS, b, 0), name=str(MIO_BREP_BREPS, b, 0),
              faces=a(MIO_BREP_BREPS, b, 0) .+ 1, edges=a(MIO_BREP_BREPS, b, 1) .+ 1,
              vertices=a(MIO_BREP_BREPS, b, 2) .+ 1, extra=str(MIO_BREP_BREPS, b, 1))
             for b in 0:n(MIO_BREP_BREPS)-1]
    (patches=patches, breps=breps, faces=faces, edges=edges, vertices=vertices,
     extra=str(MIO_BREP_FILE, 0, 0))
end

function _mdpa_data(h, section, i)
    n = Int(_check_count(ccall(_sym(:mio_mdpa_info_data_count), Int64,
                               (Ptr{Cvoid}, Int32, Int64), h, section, Int64(i)), "mdpa data"))
    out = Pair{String,Any}[]
    for e in 0:n-1
        kind = ccall(_sym(:mio_mdpa_info_data_kind), Int32, (Ptr{Cvoid}, Int32, Int64, Int64),
                     h, section, Int64(i), Int64(e))
        _check_count(kind, "mdpa data kind")
        key = _getstring((buf, len) -> ccall(_sym(:mio_mdpa_info_data_string), Int64,
                                             (Ptr{Cvoid}, Int32, Int64, Int64, Int32, Ptr{UInt8},
                                              Int64),
                                             h, section, Int64(i), Int64(e), Int32(0), buf, len))
        value = if kind == 1
            _getstring((buf, len) -> ccall(_sym(:mio_mdpa_info_data_string), Int64,
                                           (Ptr{Cvoid}, Int32, Int64, Int64, Int32, Ptr{UInt8},
                                            Int64),
                                           h, section, Int64(i), Int64(e), Int32(1), buf, len))
        else
            a = _mdpa_copy((d, t, nd, s) -> ccall(_sym(:mio_mdpa_info_data_array), Cint,
                                                  (Ptr{Cvoid}, Int32, Int64, Int64,
                                                   Ptr{Ptr{Cvoid}}, Ptr{Cint}, Ptr{Int32},
                                                   Ptr{Int64}),
                                                  h, section, Int64(i), Int64(e), d, t, nd, s))
            kind == 0 && length(a) == 1 ? a[1] : permutedims(a)
        end
        push!(out, key => value)
    end
    out
end

"""
    mdpa_info(info) -> NamedTuple

**Copy** of an MDPA side channel as plain Julia values:

* `properties` — `(id, values)` per `Begin Properties` block;
* `entity_names` — `(name, is_condition)` per cell block;
* `skipped` — what `ReadOptions(lenient=true)` skipped;
* `model_part_data` — the non-numeric `ModelPartData` entries (`key => text`);
* `tables` — `(header, values)`, `values` a `(rows, columns)` matrix;
* `geometries` — `(name, type, connectivity, ids)`, `connectivity` 1-based
  `(nodes_per_geometry, num_geometries)` like [`connectivity`](@ref);
* `mesh_blocks` — `(id, data, nodes, element_ids, condition_ids)`, `nodes`
  1-based points, the ids as the file spelled them;
* `submodelparts` — `(name, data, tables, geometry_ids, constraint_ids)`;
* `raw_blocks` — `(header, body, terminator)`, kept verbatim (`Constraints`, ...).

`data`/`values` entries are `key => value` pairs: a number, a string, or a
`(rows, columns)` matrix for an inline table.
"""
function mdpa_info(info::FormatInfo)
    h = _handle(info)
    format_name(info) == "mdpa" ||
        throw(MeshioError(MIO_ERR_INVALID_ARG, "not an mdpa side channel"))
    n(section) = _mdpa_count(h, section)
    (
        properties=[(id=_mdpa_int(h, MIO_MDPA_PROPERTIES, i, 0),
                     values=_mdpa_data(h, MIO_MDPA_PROPERTIES, i))
                    for i in 0:n(MIO_MDPA_PROPERTIES)-1],
        entity_names=[(name=_mdpa_string(h, MIO_MDPA_ENTITY_NAMES, i, 0),
                       is_condition=_mdpa_int(h, MIO_MDPA_ENTITY_NAMES, i, 0) == 1)
                      for i in 0:n(MIO_MDPA_ENTITY_NAMES)-1],
        skipped=[_mdpa_string(h, MIO_MDPA_SKIPPED, i, 0) for i in 0:n(MIO_MDPA_SKIPPED)-1],
        model_part_data=_mdpa_data(h, MIO_MDPA_MODEL_PART_DATA, 0),
        tables=[(header=_mdpa_string(h, MIO_MDPA_TABLES, i, 0),
                 values=permutedims(_mdpa_array(h, MIO_MDPA_TABLES, i, 0)))
                for i in 0:n(MIO_MDPA_TABLES)-1],
        geometries=[(name=_mdpa_string(h, MIO_MDPA_GEOMETRIES, i, 0),
                     type=_mdpa_string(h, MIO_MDPA_GEOMETRIES, i, 1),
                     connectivity=_mdpa_array(h, MIO_MDPA_GEOMETRIES, i, 0) .+ 1,
                     ids=_mdpa_array(h, MIO_MDPA_GEOMETRIES, i, 1))
                    for i in 0:n(MIO_MDPA_GEOMETRIES)-1],
        mesh_blocks=[(id=_mdpa_int(h, MIO_MDPA_MESH_BLOCKS, i, 0),
                      data=_mdpa_data(h, MIO_MDPA_MESH_BLOCKS, i),
                      nodes=_mdpa_array(h, MIO_MDPA_MESH_BLOCKS, i, 0) .+ 1,
                      element_ids=_mdpa_array(h, MIO_MDPA_MESH_BLOCKS, i, 1),
                      condition_ids=_mdpa_array(h, MIO_MDPA_MESH_BLOCKS, i, 2))
                     for i in 0:n(MIO_MDPA_MESH_BLOCKS)-1],
        submodelparts=[(name=_mdpa_string(h, MIO_MDPA_SUBMODELPARTS, i, 0),
                        data=_mdpa_data(h, MIO_MDPA_SUBMODELPARTS, i),
                         tables=_mdpa_array(h, MIO_MDPA_SUBMODELPARTS, i, 0),
                         geometry_ids=_mdpa_array(h, MIO_MDPA_SUBMODELPARTS, i, 1),
                         constraint_ids=_mdpa_array(h, MIO_MDPA_SUBMODELPARTS, i, 2))
                       for i in 0:n(MIO_MDPA_SUBMODELPARTS)-1],
        raw_blocks=[(header=_mdpa_string(h, MIO_MDPA_RAW_BLOCKS, i, 0),
                     body=_mdpa_string(h, MIO_MDPA_RAW_BLOCKS, i, 1),
                     terminator=_mdpa_string(h, MIO_MDPA_RAW_BLOCKS, i, 2))
                    for i in 0:n(MIO_MDPA_RAW_BLOCKS)-1],
    )
end

# --- the isogeometric layer (v16.43.0; doc/iga.md) -----------------------------

function _iga_opts(; subdivisions=1, elevate=0, sample_subdivisions=0, trim_samples=16,
                   edges=true, trim=true)
    three(v) = v isa Integer ? ntuple(_ -> Int32(v), 3) :
               ntuple(d -> Int32(v[min(d, length(v))]), 3)
    _CIgaOpts(three(subdivisions), three(elevate), Int32(sample_subdivisions),
              Int32(trim_samples), Int32(edges), Int32(trim), (Int32(0), Int32(0)),
              ntuple(_ -> Int64(0), 6))
end

"""
    iga_refine(info; subdivisions=1, elevate=0, sample_subdivisions=0,
               trim_samples=16, edges=true, trim=true) -> (Mesh, FormatInfo)

Refine every patch (every face surface of a B-rep) of a NURBS model exactly:
degree elevation by `elevate`, then knot insertion splitting each knot span in
`subdivisions` parts (an integer, or one per direction; an `"mfem"` model refines
every direction alike). The geometry, fields and trimming loops do not change.
`info` comes from [`read_with_info`](@ref) of an MFEM NURBS mesh, Kratos CAD JSON,
STEP or IGES. Returns the refined model sampled (`sample_subdivisions` cells per
knot span, 0 automatic) and the refined model, for [`write_with_info`](@ref). See
`doc/iga.md`.
"""
function iga_refine(info::FormatInfo; kwargs...)
    opts = _iga_opts(; kwargs...)
    refined = Ref{Ptr{Cvoid}}(C_NULL)
    ptr = ccall(_sym(:mio_iga_refine), Ptr{Cvoid},
                (Ptr{Cvoid}, Ref{_CIgaOpts}, Ptr{Ptr{Cvoid}}), _handle(info), Ref(opts), refined)
    mesh = Mesh(_check_ptr(ptr))
    mesh, FormatInfo(refined[])
end

"""
    iga_sample(info; sample_subdivisions=0, trim_samples=16, edges=true, trim=true)
        -> (Mesh, FormatInfo)

Sample a NURBS model into Lagrange cells: `sample_subdivisions` cells per knot
span and direction (0 automatic; the MFEM reader's own sampling for an `"mfem"`
model). The second value is the model with the sampling's cell-to-patch maps.
"""
function iga_sample(info::FormatInfo; kwargs...)
    opts = _iga_opts(; kwargs...)
    mapped = Ref{Ptr{Cvoid}}(C_NULL)
    ptr = ccall(_sym(:mio_iga_sample), Ptr{Cvoid},
                (Ptr{Cvoid}, Ref{_CIgaOpts}, Ptr{Ptr{Cvoid}}), _handle(info), Ref(opts), mapped)
    mesh = Mesh(_check_ptr(ptr))
    mesh, FormatInfo(mapped[])
end

# --- mesh generation from CAD (v17.13.0; doc/mesh_cad.md) ----------------------

const _MESH_CAD_FACE_STATUS_NAMES = ("ok", "empty", "open_loop", "self_intersecting",
                                     "degenerate", "projection_failed", "budget_exceeded",
                                     "unsupported", "other")

const _MESH_CAD_SMOOTH_METHODS = ("laplacian", "odt")

"""
    mesh_cad_surface(info; target_size, max_chordal_deviation=0.0,
                     min_points_per_edge=0, tolerance=0.0, record_face_ids=true,
                     edges=true, match_unreferenced_trims=true, strict=false,
                     max_points=50_000_000, min_angle=0.0, smooth_iterations=0,
                     smooth_method="laplacian", order=1)
        -> (; mesh, num_points, num_triangles, ..., faces, warnings,
              angle_histogram, corners)

Triangulate the faces of a B-rep (`info` from [`read_with_info`](@ref) of Kratos
CAD JSON, STEP or IGES) into one conforming surface mesh with edges of about
`target_size`: every edge sampled once and shared by the faces that use it, each
face a constrained Delaunay triangulation of its parameters lifted onto its
surface. A face that cannot be meshed is left out and listed in `faces` (its
`status` and `reason`), never meshed wrong; `strict=true` throws instead.
`min_angle` (degrees, under 60) refines every face until its triangles'
smallest angles on the surface reach it; the corners under 60 degrees and the
poles that keep triangles below it are in `corners` (1-based `face` and
`point`), every angle in the 18 ten-degree bins of `angle_histogram`.
`smooth_iterations` passes of guarded smoothing (`"laplacian"` or `"odt"`)
follow. `order` 2 or 3 writes curved `triangle6` or `triangle10` elements (and
`line3` or `line4` edges) whose nodes lie on the surfaces and edges; an element
whose Jacobian is not positive everywhere is kept straight-sided and counted in
`num_curved_straightened`. See `doc/mesh_cad.md`.
"""
function mesh_cad_surface(info::FormatInfo; target_size::Real, max_chordal_deviation::Real=0.0,
                          min_points_per_edge::Integer=0, tolerance::Real=0.0,
                          record_face_ids::Bool=true, edges::Bool=true,
                          match_unreferenced_trims::Bool=true, strict::Bool=false,
                          max_points::Integer=50_000_000, min_angle::Real=0.0,
                          smooth_iterations::Integer=0,
                          smooth_method::AbstractString="laplacian", order::Integer=1)
    method = findfirst(==(smooth_method), _MESH_CAD_SMOOTH_METHODS)
    method === nothing &&
        throw(ArgumentError("meshio++: mesh_cad_surface: smooth_method must be " *
                            "\"laplacian\" or \"odt\", got \"$smooth_method\""))
    opts = _CMeshCadSurfaceOpts(Cdouble(target_size), Cdouble(max_chordal_deviation),
                                Cdouble(tolerance), Int64(min_points_per_edge),
                                Int64(max_points), Int32(record_face_ids), Int32(edges),
                                Int32(match_unreferenced_trims), Int32(strict),
                                Cdouble(min_angle), Int32(smooth_iterations), Int32(method - 1),
                                Int32(order), Int32(0), ntuple(_ -> Int64(0), 3))
    report = Ref{_CMeshCadSurfaceReport}()
    ptr = ccall(_sym(:mio_mesh_cad_surface), Ptr{Cvoid},
                (Ptr{Cvoid}, Ref{_CMeshCadSurfaceOpts}, Ref{_CMeshCadSurfaceReport}),
                _handle(info), Ref(opts), report)
    mesh = Mesh(_check_ptr(ptr))
    r = report[]
    faces = map(0:r.num_faces-1) do f
        st = ccall(_sym(:mio_mesh_cad_surface_face_status), Int32, (Int64,), f)
        reason = _getstring((buf, n) -> ccall(_sym(:mio_mesh_cad_surface_face_reason), Int64,
                                              (Int64, Ptr{UInt8}, Int64), f, buf, n))
        (face=Int(f) + 1,
         status=0 <= st < length(_MESH_CAD_FACE_STATUS_NAMES) ?
                _MESH_CAD_FACE_STATUS_NAMES[st + 1] : string(st),
         reason=reason)
    end
    text = _getstring((buf, n) -> ccall(_sym(:mio_mesh_cad_surface_warnings), Int64,
                                        (Ptr{UInt8}, Int64), buf, n))
    histogram = zeros(Int64, 18)
    ccall(_sym(:mio_mesh_cad_surface_angle_histogram), Int64, (Ptr{Int64}, Int64),
          histogram, Int64(length(histogram)))
    ncorners = ccall(_sym(:mio_mesh_cad_surface_num_corners), Int64, ())
    corners = map(0:ncorners-1) do k
        c = Ref{_CMeshCadCorner}()
        _check(ccall(_sym(:mio_mesh_cad_surface_corner), Cint, (Int64, Ref{_CMeshCadCorner}),
                     k, c))
        (face=Int(c[].face) + 1, point=Int(c[].point) + 1, u=Float64(c[].u),
         v=Float64(c[].v), x=collect(Float64, c[].x), angle_deg=Float64(c[].angle_deg),
         min_forced_deg=Float64(c[].min_forced_deg), reach=Float64(c[].reach),
         pole=c[].pole != 0)
    end
    (mesh=mesh, num_points=Int(r.num_points), num_triangles=Int(r.num_triangles),
     num_edge_segments=Int(r.num_edge_segments), num_faces=Int(r.num_faces),
     num_failed=Int(r.num_failed), num_edges=Int(r.num_edges),
     num_open_edges=Int(r.num_open_edges), num_non_manifold_edges=Int(r.num_non_manifold_edges),
     num_misoriented_edges=Int(r.num_misoriented_edges),
     num_degenerate_dropped=Int(r.num_degenerate_dropped), tolerance=Float64(r.tolerance),
     min_angle_deg=Float64(r.min_angle_deg), max_angle_deg=Float64(r.max_angle_deg),
     min_edge=Float64(r.min_edge), max_edge=Float64(r.max_edge),
     max_edge_gap=Float64(r.max_edge_gap), faces=faces,
     warnings=isempty(text) ? String[] : String.(Base.split(text, '\n')),
     num_refinement_points=Int(r.num_refinement_points),
     num_faces_below_bound=Int(r.num_faces_below_bound), num_smoothed=Int(r.num_smoothed),
     min_input_angle_deg=Float64(r.min_input_angle_deg), angle_histogram=histogram,
     corners=corners, num_curved_straightened=Int(r.num_curved_straightened))
end

const MIO_IGA_POINTS = Int32(0)
const MIO_IGA_DERIVATIVES = Int32(1)
const MIO_IGA_NORMALS = Int32(2)
const MIO_IGA_INSIDE = Int32(3)
const MIO_IGA_FIELD = Int32(4)

"""
    iga_evaluate(info, patch, params; derivatives=0, trim_samples=16)
        -> (; points, fields, derivatives, normals, inside)

Evaluate patch `patch` (1-based; the face of a B-rep) of a NURBS model at the
parametric coordinates in `params`, a `(d, N)` matrix with one column per point
(`d` the parametric dimension), each inside the patch's knot range. `points` is
`(space dimension, N)`, `fields` one `(components, N)` matrix per NURBS field.
With `derivatives=1`, `derivatives` is `(space dimension, d, N)` and `normals`
`(3, N)` (surfaces in 3-D); for a B-rep face `inside` is a `Vector{Bool}` of the
trimming test. Absent parts are `nothing`.
"""
function iga_evaluate(info::FormatInfo, patch::Integer, params::AbstractMatrix{<:Real};
                      derivatives::Integer=0, trim_samples::Integer=16)
    flat = Vector{Float64}(vec(params))
    ev = ccall(_sym(:mio_iga_evaluate), Ptr{Cvoid},
               (Ptr{Cvoid}, Int64, Ptr{Float64}, Int64, Int32, Int32),
               _handle(info), Int64(patch - 1), flat, Int64(size(params, 2)),
               Int32(derivatives), Int32(trim_samples))
    _check_ptr(ev)
    try
        arr(which, index=0) = _mdpa_copy((d, t, nd, s) -> ccall(
            _sym(:mio_iga_evaluation_array), Cint,
            (Ptr{Cvoid}, Int32, Int64, Ptr{Ptr{Cvoid}}, Ptr{Cint}, Ptr{Int32}, Ptr{Int64}),
            ev, which, Int64(index), d, t, nd, s))
        opt(which) = try
            arr(which)
        catch e
            e isa MeshioError ? nothing : rethrow()
        end
        nf = Int(ccall(_sym(:mio_iga_evaluation_num_fields), Int64, (Ptr{Cvoid},), ev))
        inside = opt(MIO_IGA_INSIDE)
        (points=arr(MIO_IGA_POINTS), fields=[arr(MIO_IGA_FIELD, k) for k in 0:nf-1],
         derivatives=opt(MIO_IGA_DERIVATIVES), normals=opt(MIO_IGA_NORMALS),
         inside=inside === nothing ? nothing : Vector{Bool}(inside .!= 0))
    finally
        ccall(_sym(:mio_iga_evaluation_free), Cvoid, (Ptr{Cvoid},), ev)
    end
end
