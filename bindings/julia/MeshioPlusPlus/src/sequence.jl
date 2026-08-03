# Multi-file / transient datasets: an ordered plan over a set of files (or the
# steps inside one multi-step file), read one step at a time.
#
# See doc/sequences.md. The C-side handle stores paths, per-file step indices
# and time values, and never a mesh -- which is why `read_step` returns an
# owning `Mesh` rather than a borrow the sequence would have to cache.

"""
    Sequence(pattern; format="", times=nothing, time_from="auto")
    Sequence(paths::AbstractVector{<:AbstractString}; ..., sort=false)

An ordered, lazily-read multi-file dataset.

Built from a glob `pattern` (`*` and `?` only -- **no** `**` or `[set]`; the
directory part is taken literally) or from an explicit list of `paths`, whose
order is the caller's and is kept unless `sort=true`.

Ordering is **natural-numeric**, so `out_9.vtu` precedes `out_10.vtu` — a
lexicographic sort gets that backwards, which is why the rule is documented
rather than left to the filesystem.

Each entry carries a time value and a record of where it came from
([`time_source`](@ref)): `:explicit`, `:file`, `:filename` or `:index` (the
fallback). "The file said 0.25" and "nothing said anything, so this is
position 3" are different facts.

Nothing is read until [`read_step`](@ref), and the sequence caches nothing, so
a 500-step dataset is traversable without materialising it:

```julia
seq = Sequence("out_*.vtu")
for i in 1:length(seq)
    mesh = read_step(seq, i)
    # ... one mesh alive at a time ...
end
close(seq)
```

Like [`Mesh`](@ref), the handle is released by a finalizer; [`close`](@ref)
does it deterministically and is idempotent. The `do`-block form closes it even
if the body throws.
"""
mutable struct Sequence
    ptr::Ptr{Cvoid}

    function Sequence(ptr::Ptr{Cvoid})
        s = new(ptr)
        finalizer(s) do x
            if x.ptr != C_NULL
                ccall(_sym(:mio_sequence_free), Cvoid, (Ptr{Cvoid},), x.ptr)
                x.ptr = C_NULL
            end
        end
        s
    end
end

# Mirror of C `mio_sequence_opts`; field order and types are ABI.
struct _SequenceOpts
    format::Cstring
    times::Ptr{Cdouble}
    num_times::Int64
    time_from::Int32
    sort::Int32
    reserved::NTuple{6,Int64}
end

const _TIME_FROM = Dict("auto" => Int32(0), "file" => Int32(1),
                        "filename" => Int32(2), "index" => Int32(3))
const _TIME_SOURCE = (:explicit, :file, :filename, :index)

function _sequence_opts(format, times, time_from, sort)
    code = get(_TIME_FROM, String(time_from), nothing)
    code === nothing && throw(ArgumentError(
        "time_from must be \"auto\", \"file\", \"filename\" or \"index\", got $(time_from)"))
    tvec = times === nothing ? Cdouble[] : Cdouble[Float64(t) for t in times]
    (_SequenceOpts(Base.unsafe_convert(Cstring, Base.cconvert(Cstring, String(format))),
                   isempty(tvec) ? Ptr{Cdouble}(C_NULL) : pointer(tvec),
                   Int64(length(tvec)), code, Int32(sort), ntuple(_ -> Int64(0), 6)),
     tvec)
end

function Sequence(pattern::AbstractString; format::AbstractString="",
                  times=nothing, time_from="auto")
    fmt = Base.cconvert(Cstring, String(format))
    opts, tvec = _sequence_opts(format, times, time_from, false)
    GC.@preserve fmt tvec begin
        Sequence(_check_ptr(ccall(_sym(:mio_sequence_open_ex), Ptr{Cvoid},
                                  (Cstring, Ref{_SequenceOpts}), pattern, Ref(opts))))
    end
end

function Sequence(paths::AbstractVector{<:AbstractString}; format::AbstractString="",
                  times=nothing, time_from="auto", sort::Bool=false)
    isempty(paths) && throw(ArgumentError("the path list is empty"))
    fmt = Base.cconvert(Cstring, String(format))
    opts, tvec = _sequence_opts(format, times, time_from, sort)
    cpaths = [Base.cconvert(Cstring, String(p)) for p in paths]
    ptrs = [Base.unsafe_convert(Cstring, p) for p in cpaths]
    GC.@preserve fmt tvec cpaths begin
        Sequence(_check_ptr(ccall(_sym(:mio_sequence_open_list), Ptr{Cvoid},
                                  (Ptr{Cstring}, Int64, Ref{_SequenceOpts}),
                                  ptrs, Int64(length(ptrs)), Ref(opts))))
    end
end

function Sequence(f::Function, args...; kwargs...)
    seq = Sequence(args...; kwargs...)
    try
        f(seq)
    finally
        close(seq)
    end
end

function _handle(s::Sequence)
    s.ptr == C_NULL && throw(MeshioError(MIO_ERR_INVALID_ARG, "sequence has been closed"))
    s.ptr
end

"""
    close(seq::Sequence)

Release the sequence handle deterministically. Idempotent; the finalizer does
the same thing later if you do not.
"""
function Base.close(s::Sequence)
    if s.ptr != C_NULL
        ccall(_sym(:mio_sequence_free), Cvoid, (Ptr{Cvoid},), s.ptr)
        s.ptr = C_NULL
    end
    nothing
end

"""
    length(seq::Sequence)

The number of steps.
"""
Base.length(s::Sequence) =
    Int(_check_count(ccall(_sym(:mio_sequence_count), Int64, (Ptr{Cvoid},), _handle(s)),
                     "mio_sequence_count"))

"""
    path(seq, i)

Entry `i`'s file path (1-based).
"""
function path(s::Sequence, i::Integer)
    handle = _handle(s)
    _getstring() do buf, len
        ccall(_sym(:mio_sequence_path), Int64, (Ptr{Cvoid}, Int64, Ptr{UInt8}, Int64),
              handle, Int64(i - 1), buf, len)
    end
end

"""
    step(seq, i)

Entry `i`'s step index **within its own file** — 0 for a single-step file, so a
directory of `.vtu`s and one multi-step `.xdmf` have the same shape.
"""
step(s::Sequence, i::Integer) =
    Int(_check_count(ccall(_sym(:mio_sequence_step), Int64, (Ptr{Cvoid}, Int64),
                           _handle(s), Int64(i - 1)), "mio_sequence_step"))

"""
    time(seq, i)

Entry `i`'s time value.
"""
function time(s::Sequence, i::Integer)
    out = Ref{Cdouble}(0.0)
    _check(ccall(_sym(:mio_sequence_time), Cint, (Ptr{Cvoid}, Int64, Ref{Cdouble}),
                 _handle(s), Int64(i - 1), out))
    out[]
end

"""
    time_source(seq, i)

Where entry `i`'s time came from: `:explicit`, `:file`, `:filename` or
`:index`.
"""
function time_source(s::Sequence, i::Integer)
    code = ccall(_sym(:mio_sequence_time_source), Int32, (Ptr{Cvoid}, Int64),
                 _handle(s), Int64(i - 1))
    _check_count(code, "mio_sequence_time_source")
    _TIME_SOURCE[code + 1]
end

"""
    read_step(seq, i)

Read entry `i` (1-based). The returned [`Mesh`](@ref) is **owned** by the
caller and independent of the sequence, which caches nothing — that is what
keeps a long dataset traversable.
"""
read_step(s::Sequence, i::Integer) =
    Mesh(_check_ptr(ccall(_sym(:mio_sequence_read), Ptr{Cvoid}, (Ptr{Cvoid}, Int64),
                          _handle(s), Int64(i - 1))))

# Iteration: `for mesh in seq` yields one owned mesh per step, lazily.
Base.iterate(s::Sequence, i::Int=1) = i > length(s) ? nothing : (read_step(s, i), i + 1)
Base.eltype(::Type{Sequence}) = Mesh
Base.IteratorSize(::Type{Sequence}) = Base.HasLength()

"""
    to_timeseries(seq, out_path; format="")

Fan-in: write every step into one multi-step file. Streams — one mesh alive at
a time. A format that cannot hold a series throws naming itself and pointing at
`{step}`, never a silent truncation to the first step.
"""
function to_timeseries(s::Sequence, out_path::AbstractString; format::AbstractString="")
    _check(ccall(_sym(:mio_sequence_to_timeseries), Cint, (Ptr{Cvoid}, Cstring, Cstring),
                 _handle(s), out_path, format))
    nothing
end

"""
    timeseries_to_sequence(in_path, out_pattern; in_format="", out_format="")

Fan-out: write each step of the multi-step file `in_path` to `out_pattern`,
which must contain `{step}` or `{index}`. Streams likewise.
"""
function timeseries_to_sequence(in_path::AbstractString, out_pattern::AbstractString;
                                in_format::AbstractString="", out_format::AbstractString="")
    _check(ccall(_sym(:mio_timeseries_to_sequence), Cint,
                 (Cstring, Cstring, Cstring, Cstring),
                 in_path, in_format, out_pattern, out_format))
    nothing
end

"""
    run_sequence_file(settings_path)

Run a sequence settings document: the pipeline schema plus `Mode`,
`Input.Pattern`, `Input.Paths`, `Input.Times` and `Input.TimeFrom`. A document
using none of those behaves exactly as [`run_pipeline_file`](@ref).
"""
function run_sequence_file(settings_path::AbstractString)
    _check(ccall(_sym(:mio_sequence_pipeline_run_file), Cint, (Cstring,), settings_path))
    nothing
end

"""
    run_sequence_json(json_text)

[`run_sequence_file`](@ref) over JSON text.
"""
function run_sequence_json(json_text::AbstractString)
    _check(ccall(_sym(:mio_sequence_pipeline_run_json), Cint, (Cstring,), json_text))
    nothing
end
