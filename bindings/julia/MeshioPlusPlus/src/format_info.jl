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

"""
    FormatInfo

A format's side channel, returned by [`read_with_info`](@ref). Released by a
finalizer; [`close`](@ref) does it deterministically and is idempotent. Inspect
an MDPA one with [`mdpa_info`](@ref), or a Gmsh one with [`gmsh_info`](@ref).
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
(currently `"mdpa"` and `"gmsh"`) the second value is a [`FormatInfo`](@ref); for every other
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

The format the side channel belongs to (`"mdpa"` or `"gmsh"`).
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
