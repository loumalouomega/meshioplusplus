# Notebook-only display helpers for the meshio++ Julia examples.
#
# This binding rides on the flat C API (see doc/julia.md), and CLAUDE.md
# documents that data-driven SVG colouring is a flat-ABI gap: registry.cpp's
# (path, mesh) writer lambdas can't carry per-call parameters, so `mio.write`
# always emits the fixed default styling -- no color_by, no colorbar, no
# camera control (unlike the C++ notebooks, which call `write_svg` directly
# and bypass the registry). So where the C++ notebooks colour a render by a
# quality metric or a field, these notebooks show the same information as a
# small hand-rolled SVG chart instead -- honest about the gap rather than
# working around it.
#
# Nothing here is part of the package API; it exists only to give these
# notebooks something to `display()`.

using MeshioPlusPlus
import MeshioPlusPlus as mio

# -- rich display plumbing ----------------------------------------------------

"""A raw SVG string, displayed inline via Jupyter's `image/svg+xml` MIME type."""
struct SvgImage
    svg::String
end

Base.show(io::IO, ::MIME"image/svg+xml", x::SvgImage) = print(io, x.svg)
Base.showable(::MIME"image/svg+xml", ::SvgImage) = true

let _n = Ref(0)
    global next_temp_path
    """A fresh scratch file path with the given suffix, for one-shot round trips."""
    function next_temp_path(suffix::AbstractString)
        _n[] += 1
        joinpath(tempdir(), "meshioplusplus_julia_nb_$(_n[])$suffix")
    end
end

# -- mesh renders, via meshio++'s own SVG writer ------------------------------

"""
    render(mesh) -> SvgImage

Write `mesh` through the C API's SVG writer (fixed default camera and
styling -- see the module note above) and return it as a displayable
[`SvgImage`](@ref).
"""
function render(m::Mesh)
    path = next_temp_path(".svg")
    mio.write(m, path; format="svg")
    SvgImage(read(path, String))
end

# -- synthetic fixtures --------------------------------------------------------

"""
    hex_block(n; spacing=1.0) -> Mesh

An `n`x`n`x`n` hexahedron block on the lattice `[0, n*spacing]^3` -- the
small, easy-to-read fixture several operation demos build on
(`convert_cells`, `refine`, `smooth`, `interpolate`), mirroring the C++
notebooks' `hex_block` and the Python notebooks' `np.meshgrid`-based one.
"""
function hex_block(n::Integer; spacing::Real=1.0)
    pid(i, j, k) = (i * (n + 1) + j) * (n + 1) + k + 1  # 1-based flat index
    npts = (n + 1)^3
    pts = Matrix{Float64}(undef, 3, npts)
    for i in 0:n, j in 0:n, k in 0:n
        id = pid(i, j, k)
        pts[1, id] = i * spacing
        pts[2, id] = j * spacing
        pts[3, id] = k * spacing
    end

    ncells = n^3
    conn = Matrix{Int64}(undef, 8, ncells)
    cell = 0
    for i in 0:(n-1), j in 0:(n-1), k in 0:(n-1)
        cell += 1
        conn[:, cell] = [pid(i, j, k), pid(i + 1, j, k), pid(i + 1, j + 1, k), pid(i, j + 1, k),
                        pid(i, j, k + 1), pid(i + 1, j, k + 1), pid(i + 1, j + 1, k + 1),
                        pid(i, j + 1, k + 1)]
    end

    m = Mesh()
    set_points!(m, pts)
    add_cell_block!(m, "hexahedron", conn)
    m
end

# -- small hand-rolled charts (no plotting library; see the module note) -----

_svg_escape(s::AbstractString) = replace(replace(replace(s, "&" => "&amp;"), "<" => "&lt;"), ">" => "&gt;")

"""
    bar_chart(labels, values, title; color="#3b82f6") -> SvgImage

Horizontal bar chart -- e.g. file sizes per format, cell counts per type.
Mirrors the C++ notebooks' `bar_chart` helper (`mio_notebook.hpp`).
"""
function bar_chart(labels::Vector{<:AbstractString}, values::Vector{<:Real}, title::AbstractString;
                   color::AbstractString="#3b82f6")
    n = length(labels)
    barH, gap, left, right, top, width = 30.0, 12.0, 220.0, 60.0, 50.0, 760.0
    maxv = maximum(values)
    plotw = width - left - right
    height = top + n * (barH + gap) + 20.0

    io = IOBuffer()
    print(io, "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 $width $height\" ",
          "font-family=\"sans-serif\" font-size=\"13\">")
    print(io, "<text x=\"$(width/2)\" y=\"24\" text-anchor=\"middle\" font-size=\"16\">",
          _svg_escape(title), "</text>")
    for i in 1:n
        y = top + (i - 1) * (barH + gap)
        w = maxv > 0 ? (values[i] / maxv) * plotw : 0.0
        print(io, "<text x=\"$(left-10)\" y=\"$(y+barH*0.65)\" text-anchor=\"end\">",
              _svg_escape(labels[i]), "</text>")
        print(io, "<rect x=\"$left\" y=\"$y\" width=\"$w\" height=\"$barH\" fill=\"$color\" />")
        print(io, "<text x=\"$(left+w+8)\" y=\"$(y+barH*0.65)\">", round(values[i]; digits=4), "</text>")
    end
    print(io, "</svg>")
    SvgImage(String(take!(io)))
end

"""
    histogram_chart(edges, counts, title, xlabel; color="#3b82f6") -> SvgImage

Vertical histogram -- `length(edges) == length(counts) + 1`. Mirrors the
C++ notebooks' `histogram_chart` helper.
"""
function histogram_chart(edges::Vector{<:Real}, counts::Vector{<:Integer}, title::AbstractString,
                         xlabel::AbstractString; color::AbstractString="#3b82f6")
    n = length(counts)
    width, height, left, right, top, bottom = 760.0, 380.0, 60.0, 20.0, 50.0, 60.0
    plotw, ploth = width - left - right, height - top - bottom
    maxc = maximum(counts)
    barw = plotw / n

    io = IOBuffer()
    print(io, "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 $width $height\" ",
          "font-family=\"sans-serif\" font-size=\"12\">")
    print(io, "<text x=\"$(width/2)\" y=\"24\" text-anchor=\"middle\" font-size=\"16\">",
          _svg_escape(title), "</text>")
    print(io, "<line x1=\"$left\" y1=\"$(top+ploth)\" x2=\"$(left+plotw)\" y2=\"$(top+ploth)\" ",
          "stroke=\"black\" />")
    for i in 1:n
        h = maxc > 0 ? (counts[i] / maxc) * ploth : 0.0
        x = left + (i - 1) * barw
        print(io, "<rect x=\"$(x+barw*0.05)\" y=\"$(top+ploth-h)\" width=\"$(barw*0.9)\" height=\"$h\" ",
              "fill=\"$color\" />")
    end
    print(io, "<text x=\"$left\" y=\"$(top+ploth+20)\">", round(first(edges); digits=3), "</text>")
    print(io, "<text x=\"$(left+plotw)\" y=\"$(top+ploth+20)\" text-anchor=\"end\">",
          round(last(edges); digits=3), "</text>")
    print(io, "<text x=\"$(width/2)\" y=\"$(height-10)\" text-anchor=\"middle\">", _svg_escape(xlabel), "</text>")
    print(io, "</svg>")
    SvgImage(String(take!(io)))
end
