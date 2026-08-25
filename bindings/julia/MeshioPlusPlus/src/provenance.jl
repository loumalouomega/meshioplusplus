# Provenance (v10.16.0 write side, v10.17.0 read side) -- see doc/provenance.md.
#
# Rides the flat C ABI (`mio_provenance_*`) exactly as every other operation
# here does. The one shaping choice: the C ABI is a begin/end pair (a scope
# cannot cross an ABI boundary as RAII), but exposing that raw would make a
# thrown exception inside the block leak the scope. `provenance` is therefore
# a `do`-block form with the end call in a `finally`, matching how `Sequence`
# and `XdmfSeries` already give C handles a Julia-shaped lifetime.

"""
    ProvenanceMode

How much provenance a writer should render: `PROVENANCE_OFF` (the default
everywhere -- just the one-line credit), `PROVENANCE_BEST_EFFORT` (the full
block where the target format's header slot allows it, degrading silently
otherwise) or `PROVENANCE_REQUIRED` (an error for a format with no slot at
all).
"""
const PROVENANCE_OFF = 0
const PROVENANCE_BEST_EFFORT = 1
const PROVENANCE_REQUIRED = 2

"""
    provenance(f, mode = PROVENANCE_BEST_EFFORT)

Run `f` with a provenance scope open, so any `write` inside it records the
richer block. The scope is closed even if `f` throws.

```julia
provenance() do
    provenance_source("bracket.msh", "gmsh")
    provenance_note("regions-dropped", "Side regions have no VTU equivalent")
    MeshioPlusPlus.write(mesh, "out.vtu")
end
```
"""
function provenance(f::Function, mode::Integer = PROVENANCE_BEST_EFFORT)
    _check(ccall(_sym(:mio_provenance_scope_begin), Cint, (Cint,), Cint(mode)))
    try
        return f()
    finally
        ccall(_sym(:mio_provenance_scope_end), Cvoid, ())
    end
end

"""
    provenance_note(category, detail)

Record one conversion assumption against the active scope. A no-op outside
one, so it is safe to call unconditionally.
"""
function provenance_note(category::AbstractString, detail::AbstractString)
    ccall(_sym(:mio_provenance_note), Cvoid, (Cstring, Cstring), category, detail)
    return nothing
end

"""
    provenance_source(path, format)

Record where the mesh came from. A no-op outside a scope.
"""
function provenance_source(path::AbstractString, format::AbstractString)
    ccall(_sym(:mio_provenance_set_source), Cvoid, (Cstring, Cstring), path, format)
    return nothing
end

"""
    provenance_target(format; encoding = "", codec = "", float_format = "")

Record what was actually written (as opposed to what was requested). A no-op
outside a scope.
"""
function provenance_target(format::AbstractString; encoding::AbstractString = "",
                           codec::AbstractString = "", float_format::AbstractString = "")
    ccall(_sym(:mio_provenance_set_target), Cvoid, (Cstring, Cstring, Cstring, Cstring),
          format, encoding, codec, float_format)
    return nothing
end
