# Software rendering (v16.33.0, field rendering v16.34.0): a deterministic
# software rasterizer and its terminal, HTML and PNG encodings. See doc/tui.md.

"""
    Frame

A rendered image, copied out of the library: `width` and `height` in pixels,
`rgba::Array{UInt8,3}` of size `(4, width, height)` (straight alpha, row 1 on
top, so `rgba[:, x, y]` is pixel `(x, y)`), `cell_ids::Matrix{Int64}` of size
`(width, height)` (the input cell drawn at each pixel, global block-major, -1 for
none), `range` (the mapped `(vmin, vmax)`, or `nothing` when no field is
mapped) and `notes` (the colour range, ticks and keys, as text).
"""
struct Frame
    width::Int
    height::Int
    rgba::Array{UInt8,3}
    cell_ids::Matrix{Int64}
    range::Union{Nothing,Tuple{Float64,Float64}}
    notes::Vector{String}
end

const _RENDER_SHADING = Dict("none" => 0, "flat" => 1, "smooth" => 2)
const _RENDER_EDGES = Dict("none" => 0, "all" => 1, "feature" => 2)
const _RENDER_SCALE = Dict("linear" => 0, "log" => 1, "symlog" => 2)
const _RENDER_DIAGNOSTIC = Dict("none" => 0, "quality" => 1, "inverted" => 2,
                                "degenerate" => 3, "orientation" => 4, "free_edges" => 5,
                                "edge_length" => 6)
const _TEXT_ENCODING = Dict("halfblock" => 0, "quadrant" => 1, "sextant" => 2, "braille" => 3,
                            "ascii" => 4, "kitty" => 5, "iterm2" => 6, "sixel" => 7)
const _TEXT_COLOR_DEPTH = Dict("truecolor" => 0, "256" => 1, "16" => 2, "mono" => 3)
const _TEXT_FORMAT = Dict("ansi" => 0, "plain" => 1, "html" => 2)

function _enum_value(table::Dict, what::AbstractString, v)
    v isa Integer && return Int32(v)
    key = String(v)
    haskey(table, key) && return Int32(table[key])
    error("meshio++: render: $what must be one of " * join(sort!(collect(keys(table))), ", ") *
          ", not '$key'")
end

"""`0xRRGGBBAA` from four channels in 0..255 (alpha defaults to opaque)."""
rgba(r::Integer, g::Integer, b::Integer, a::Integer=255) =
    (UInt32(r) << 24) | (UInt32(g) << 16) | (UInt32(b) << 8) | UInt32(a)

_color(c::Integer) = UInt32(c)
_color(c::Union{Tuple,AbstractVector}) = rgba(c...)

# Set one field of an isbits option struct that lives in a Ref.
function _setopt!(ref::Ref{T}, name::Symbol, v) where {T}
    i = Base.fieldindex(T, name)
    F = fieldtype(T, i)
    GC.@preserve ref begin
        p = Base.unsafe_convert(Ptr{T}, ref)
        unsafe_store!(Ptr{F}(p + fieldoffset(T, i)), Base.convert(F, v))
    end
    nothing
end

const _RENDER_STRINGS = (:view, :color_by, :cmap, :reduce, :expr, :vectors, :warp, :quality_metric)
const _RENDER_FLOATS = (:pixel_aspect, :azimuth, :elevation, :roll, :fov_deg, :zoom, :pan_x,
                        :pan_y, :ambient, :light_x, :light_y, :light_z, :split_angle,
                        :feature_angle, :point_radius, :scale_threshold, :vector_length,
                        :warp_scale)
const _RENDER_INTS = (:width, :height, :supersample, :component, :isolines, :vector_count)
const _RENDER_BOOLS = (:perspective, :two_sided, :colorbar, :axes, :scale_bar, :symmetric,
                       :categorical, :color_regions, :category_edges, :warp_outline)
const _RENDER_COLORS = (:edge_color, :fill_color, :line_color, :background, :nan_color,
                        :iso_color, :vector_color, :outline_color)

"""
    _with_render_opts(f; kwargs...)

Initialize a `mio_render_opts` in C, set the keyword arguments the caller named,
and call `f(Ref{_CRenderOpts})` with every borrowed string kept alive.
"""
function _with_render_opts(f; kwargs...)
    ref = Ref{_CRenderOpts}()
    GC.@preserve ref begin
        ccall(_sym(:mio_render_opts_init), Cvoid, (Ptr{_CRenderOpts},),
              Base.unsafe_convert(Ptr{_CRenderOpts}, ref))
    end
    roots = Any[]
    levels = Float64[]
    for (k, v) in kwargs
        v === nothing && continue
        if k in _RENDER_STRINGS
            bytes = Vector{UInt8}(codeunits(String(v) * "\0"))
            push!(roots, bytes)
            _setopt!(ref, k, pointer(bytes))
        elseif k in _RENDER_FLOATS
            _setopt!(ref, k, Float64(v))
        elseif k in _RENDER_INTS
            _setopt!(ref, k, Int32(v))
        elseif k in _RENDER_BOOLS
            _setopt!(ref, k, Int32(Bool(v)))
        elseif k in _RENDER_COLORS
            _setopt!(ref, k, _color(v))
        elseif k === :shading
            _setopt!(ref, k, _enum_value(_RENDER_SHADING, "shading", v))
        elseif k === :edges
            _setopt!(ref, k, _enum_value(_RENDER_EDGES, "edges", v))
        elseif k === :scale
            _setopt!(ref, k, _enum_value(_RENDER_SCALE, "scale", v))
        elseif k === :diagnostic
            _setopt!(ref, k, _enum_value(_RENDER_DIAGNOSTIC, "diagnostic", v))
        elseif k === :projection
            v in ("orthographic", "perspective") ||
                error("meshio++: render: projection must be \"orthographic\" or \"perspective\"")
            _setopt!(ref, :perspective, Int32(v == "perspective"))
        elseif k === :vmin
            _setopt!(ref, :vmin, Float64(v)); _setopt!(ref, :has_vmin, Int32(1))
        elseif k === :vmax
            _setopt!(ref, :vmax, Float64(v)); _setopt!(ref, :has_vmax, Int32(1))
        elseif k === :clip
            length(v) == 2 || error("meshio++: render: clip must be (low, high) percentiles")
            if v[1] !== nothing
                _setopt!(ref, :clip_low, Float64(v[1])); _setopt!(ref, :has_clip_low, Int32(1))
            end
            if v[2] !== nothing
                _setopt!(ref, :clip_high, Float64(v[2])); _setopt!(ref, :has_clip_high, Int32(1))
            end
        elseif k === :pan
            length(v) == 2 || error("meshio++: render: pan must be (x, y)")
            _setopt!(ref, :pan_x, Float64(v[1])); _setopt!(ref, :pan_y, Float64(v[2]))
        elseif k === :light_dir
            length(v) == 3 || error("meshio++: render: light_dir must be (x, y, z)")
            _setopt!(ref, :light_x, Float64(v[1])); _setopt!(ref, :light_y, Float64(v[2]))
            _setopt!(ref, :light_z, Float64(v[3]))
        elseif k === :iso_levels
            levels = Float64.(collect(v))
            _setopt!(ref, :iso_levels, pointer(levels))
            _setopt!(ref, :num_iso_levels, length(levels))
        else
            error("meshio++: render: unknown option '$k'")
        end
    end
    GC.@preserve roots levels f(ref)
end

function _text_opts(; encoding="halfblock", color_depth="truecolor", format="ansi", cols::Integer=80,
                    rows::Integer=24, cell_aspect::Real=2.0, cell_pixels=(8, 16),
                    tmux::Bool=false, notes::Bool=true)
    _CTextOpts(Float64(cell_aspect), _enum_value(_TEXT_ENCODING, "encoding", encoding),
               _enum_value(_TEXT_COLOR_DEPTH, "color_depth", color_depth),
               _enum_value(_TEXT_FORMAT, "format", format), Int32(cols), Int32(rows),
               Int32(cell_pixels[1]), Int32(cell_pixels[2]), Int32(tmux), Int32(notes), Int32(0),
               ntuple(_ -> Int64(0), 6))
end

"""
    render(m; width=320, height=240, kwargs...) -> Frame

Render a mesh into an RGBA frame (volume cells are drawn through their boundary
skin). The keyword arguments are those of the C `mio_render_opts`, by name:
`supersample`, `view` (`"iso"`, `"+x"` ... `"-z"`), `azimuth`, `elevation`, `roll`,
`projection`, `fov_deg`, `zoom`, `pan`, `shading` (`"none"`, `"flat"`, `"smooth"`),
`two_sided`, `ambient`, `light_dir`, `split_angle`, `edges` (`"none"`, `"all"`,
`"feature"`), `feature_angle`, the colours `edge_color`, `fill_color`, `line_color`,
`background`, `nan_color`, `iso_color`, `vector_color`, `outline_color` (a
`0xRRGGBBAA` integer or an `(r, g, b[, a])` tuple, see [`rgba`](@ref)),
`point_radius`; colouring by `color_by`, `component`, `cmap`, `vmin`, `vmax`,
`colorbar`, `axes`, `scale_bar`; and the field options `expr`, `reduce`, `clip`,
`symmetric`, `scale`, `scale_threshold`, `categorical`, `color_regions`,
`category_edges`, `isolines`, `iso_levels`, `vectors`, `vector_count`,
`vector_length`, `warp`, `warp_scale`, `warp_outline`, `diagnostic` and
`quality_metric`. See `doc/tui.md`.
"""
function render(m::Mesh; kwargs...)
    ptr = _with_render_opts(; kwargs...) do ref
        ccall(_sym(:mio_render), Ptr{Cvoid}, (Ptr{Cvoid}, Ref{_CRenderOpts}), _handle(m), ref)
    end
    h = _check_ptr(ptr)
    try
        w = Int(ccall(_sym(:mio_frame_width), Int32, (Ptr{Cvoid},), h))
        ht = Int(ccall(_sym(:mio_frame_height), Int32, (Ptr{Cvoid},), h))
        rgba_view = unsafe_wrap(Array, ccall(_sym(:mio_frame_rgba), Ptr{UInt8}, (Ptr{Cvoid},), h),
                                (4, w, ht))
        ids_view = unsafe_wrap(Array, ccall(_sym(:mio_frame_cell_ids), Ptr{Int64}, (Ptr{Cvoid},), h),
                               (w, ht))
        vmin = Ref(0.0)
        vmax = Ref(0.0)
        colored = ccall(_sym(:mio_frame_range), Int32, (Ptr{Cvoid}, Ptr{Float64}, Ptr{Float64}),
                        h, vmin, vmax)
        colored < 0 && _throw_status(MIO_ERR_INVALID_ARG)
        nnotes = Int(ccall(_sym(:mio_frame_num_notes), Int64, (Ptr{Cvoid},), h))
        notes = String[_getstring((buf, len) -> ccall(_sym(:mio_frame_note), Int64,
                                    (Ptr{Cvoid}, Int64, Ptr{UInt8}, Int64), h, i - 1, buf, len))
                       for i in 1:nnotes]
        Frame(w, ht, copy(rgba_view), copy(ids_view), colored == 1 ? (vmin[], vmax[]) : nothing,
              notes)
    finally
        ccall(_sym(:mio_frame_free), Cvoid, (Ptr{Cvoid},), h)
    end
end

"""
    render_text(m; cols=80, rows=24, encoding="halfblock", color_depth="truecolor",
                format="ansi", cell_aspect=2.0, kwargs...) -> String

Render a mesh as text sized to `cols` x `rows` terminal cells. `encoding` is
`"halfblock"`, `"quadrant"`, `"sextant"`, `"braille"`, `"ascii"`, or a graphics
protocol (`"kitty"`, `"iterm2"`, `"sixel"`, with `tmux=true` to wrap it for tmux
passthrough); `color_depth` is `"truecolor"`, `"256"`, `"16"` or `"mono"`;
`format` is `"ansi"`, `"plain"` or `"html"`. The other keyword arguments are
those of [`render`](@ref). Inside a notebook, `display("text/html", render_text(m;
format="html"))` shows the terminal form.
"""
function render_text(m::Mesh; cols::Integer=80, rows::Integer=24, encoding="halfblock",
                     color_depth="truecolor", format="ansi", cell_aspect::Real=2.0,
                     cell_pixels=(8, 16), tmux::Bool=false, notes::Bool=true, kwargs...)
    topts = _text_opts(; encoding, color_depth, format, cols, rows, cell_aspect, cell_pixels,
                       tmux, notes)
    _with_render_opts(; kwargs...) do ref
        _getstring((buf, len) -> ccall(_sym(:mio_render_text), Int64,
                       (Ptr{Cvoid}, Ref{_CRenderOpts}, Ref{_CTextOpts}, Ptr{UInt8}, Int64),
                       _handle(m), ref, Ref(topts), buf, len))
    end
end

"""
    render_png(m; compress=0, kwargs...) -> Vector{UInt8}

Render a mesh and encode it as an RGBA PNG. `compress = 0` stores the pixels in
uncompressed deflate blocks (no zlib, the same bytes everywhere); 1-9 compress
through zlib where the build has it. The keyword arguments are those of
[`render`](@ref).
"""
function render_png(m::Mesh; compress::Integer=0, kwargs...)
    ptr = _with_render_opts(; kwargs...) do ref
        ccall(_sym(:mio_render), Ptr{Cvoid}, (Ptr{Cvoid}, Ref{_CRenderOpts}), _handle(m), ref)
    end
    h = _check_ptr(ptr)
    try
        n = ccall(_sym(:mio_frame_png), Int64, (Ptr{Cvoid}, Int32, Ptr{UInt8}, Int64),
                  h, Int32(compress), Ptr{UInt8}(C_NULL), 0)
        n < 0 && _throw_status(MIO_ERR_INVALID_ARG)
        buf = Vector{UInt8}(undef, n)
        got = GC.@preserve buf ccall(_sym(:mio_frame_png), Int64,
                                     (Ptr{Cvoid}, Int32, Ptr{UInt8}, Int64),
                                     h, Int32(compress), pointer(buf), n)
        got < 0 && _throw_status(MIO_ERR_INVALID_ARG)
        buf
    finally
        ccall(_sym(:mio_frame_free), Cvoid, (Ptr{Cvoid},), h)
    end
end

"""
    write_snapshot(path, m; png_compress=0, cast_frames=36, cast_fps=12.0,
                   cast_degrees=360.0, cols=80, rows=24, encoding="halfblock",
                   color_depth="truecolor", kwargs...)

Render a mesh to a file chosen by extension: `.png`, `.txt`, `.ansi`, `.html` or
`.cast` (an asciinema orbit). Text forms are sized `cols` x `rows`; the other
keyword arguments are those of [`render`](@ref).
"""
function write_snapshot(path::AbstractString, m::Mesh; png_compress::Integer=0,
                        cast_frames::Integer=36, cast_fps::Real=12.0, cast_degrees::Real=360.0,
                        cols::Integer=80, rows::Integer=24, encoding="halfblock",
                        color_depth="truecolor", format="ansi", cell_aspect::Real=2.0,
                        cell_pixels=(8, 16), tmux::Bool=false, notes::Bool=true, kwargs...)
    topts = _text_opts(; encoding, color_depth, format, cols, rows, cell_aspect, cell_pixels,
                       tmux, notes)
    sopts = _CSnapshotOpts(Float64(cast_fps), Float64(cast_degrees), Int32(png_compress),
                           Int32(cast_frames), ntuple(_ -> Int64(0), 6))
    _with_render_opts(; kwargs...) do ref
        _check(ccall(_sym(:mio_write_snapshot), Cint,
                     (Cstring, Ptr{Cvoid}, Ref{_CRenderOpts}, Ref{_CTextOpts}, Ref{_CSnapshotOpts}),
                     String(path), _handle(m), ref, Ref(topts), Ref(sopts)))
    end
    nothing
end

"""
    detect_color_depth(; no_color=nothing, color_term=nothing, term=nothing) -> String

The colour depth a terminal advertises (`"truecolor"`, `"256"`, `"16"` or
`"mono"`) from the values of `NO_COLOR`, `COLORTERM` and `TERM`.
"""
function detect_color_depth(; no_color=nothing, color_term=nothing, term=nothing)
    a(x) = x === nothing ? Cstring(C_NULL) : x
    n = ccall(_sym(:mio_detect_color_depth), Int32, (Cstring, Cstring, Cstring),
              a(no_color), a(color_term), a(term))
    ("truecolor", "256", "16", "mono")[n + 1]
end
