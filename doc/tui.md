# Terminal rendering

meshio++ can draw a mesh where no window, browser or GPU is available: an SSH session on a cluster login node, a container, a CI log, a terminal with nothing else installed. A software rasterizer turns the mesh into pixels, and the pixels become Unicode block characters with terminal colours, an HTML page, an asciinema recording or a PNG.

```sh
meshioplusplus snapshot part.vtu -                          # draw it in this terminal
meshioplusplus snapshot part.vtu part.png --edges feature   # a PNG, no display needed
meshioplusplus snapshot result.vtu - --color-by temperature --colorbar --cmap turbo
meshioplusplus tui part.vtu                                 # orbit, zoom and pan it right here
```

```python
import meshioplusplus as mio

mesh = mio.read("part.vtu")
print(mio.render_text(mesh, color_by="temperature"))      # ANSI text for this terminal
image = mio.render_image(mesh, 800, 600, supersample=2)   # an (H, W, 4) uint8 array
mio.snapshot(mesh, "part.png", edges="all")               # .png, .txt, .ansi, .html or .cast
mio.tui(mesh)                                             # the interactive viewer
```

It is a preview, not a replacement for the [interactive viewers](./viewer.md): no GPU lighting, transparency or picking. What it offers instead is that it runs anywhere the core runs — the release binaries included, since it needs no OpenGL and vendors nothing — and that its output is reproducible to the byte.

## What is drawn

The rasterizer draws what the [SVG writer](./formats/svg.md) draws, but resolves visibility per pixel with a depth buffer instead of sorting whole faces, so faces that intersect or overlap come out right:

- 2-D cells as they are (triangles, quads, polygons, and the corners of their higher-order variants);
- volume cells, polyhedra included, through their boundary skin ([`extract_skin`](./extract_skin.md)); cell data reaches the skin through `"surface:parent_cell"`, as in the SVG writer;
- `line` cells as lines, and `vertex` cells — or every point of a mesh with no cells — as discs.

A volume mesh that also carries its own boundary faces draws both; where they coincide, the skin wins the depth tie. A quad is split on its shorter diagonal and a polygon is fanned from its first corner, so a warped quad shades the same way every time.

## Outputs

| Form | How | Notes |
| --- | --- | --- |
| this terminal | `snapshot IN -`, `render_text()` | sized to the terminal; ANSI colour escapes |
| `.ansi` | `snapshot IN out.ansi` | the same text in a file, for `cat` |
| `.txt` | `snapshot IN out.txt` | glyphs only, no escape sequence |
| `.html` | `snapshot IN out.html` | a self-contained page holding a `<pre>` of coloured spans |
| `.cast` | `snapshot IN out.cast` | an [asciicast v2](https://docs.asciinema.org/manual/asciicast/v2/) orbit around the mesh, for `asciinema play` |
| `.png` | `snapshot IN out.png`, `encode_png()` | RGBA, with no dependency |
| an array | `render_image()` | `(H, W, 4)` uint8, plus an optional id buffer |

These are renderings, not mesh files: they are not in the [format registry](./formats.md), `write()` does not produce them, and the conformance matrix does not apply.

The PNG writer stores its image data in uncompressed deflate blocks with its own CRC-32 and Adler-32, so it needs no zlib and writes the same bytes on every platform. `--png-compress 1..9` (`png_compress=` in Python) compresses through zlib instead, where the build has it; the result is smaller but can differ between zlib versions.

A `.cast` recording sweeps the azimuth through `--cast-degrees` (360 by default) in `--cast-frames` frames (36), stamped at `--cast-fps` (12) frames per second. The timestamps come from the frame index, never the clock, so the file is reproducible.

## Cell encodings

A terminal cell holds one character in one foreground and one background colour. Each encoding splits a cell into a grid of pixels and, for every cell, picks the split of its pixels into two colour groups with the least squared error, ties broken by a fixed order, then the character that draws that split:

| `--encoding` | pixels per cell | characters | best for |
| --- | --- | --- | --- |
| `halfblock` (default) | 1 x 2 | `▀ ▄ █` | colour fields: square pixels, two colours per cell |
| `quadrant` | 2 x 2 | `▘ ▝ ▖ ▗ ▚ ▞ ▙ ▟ …` | more shape detail |
| `sextant` | 2 x 3 | U+1FB00–U+1FB3B (Unicode 13) | finer shape; needs a recent font |
| `braille` | 2 x 4 | U+2800–U+28FF | line work, wireframes |
| `ascii` | 1 x 1 | ` .:-=+*#%@` by luminance | anything that cannot print Unicode |

Transparent pixels (the default background) leave the terminal's own background showing: a cell with nothing drawn is a space with no colour at all.

`--color-depth` picks the colours: `truecolor` (24-bit), `256` (the xterm palette, nearest colour, no dithering), `16` (the ANSI colours) or `mono` (no colour: a pixel is drawn when it is bright). `auto`, the default when drawing to the terminal, reads the environment: `NO_COLOR` set to anything gives `mono` (see [no-color.org](https://no-color.org)); `COLORTERM=truecolor` or `24bit` gives `truecolor`; a `TERM` containing `256color` gives `256`; `TERM=dumb` gives `mono`; anything else gives `16`. On Windows, Windows Terminal (`WT_SESSION`) counts as 24-bit. Files default to `truecolor`.

Cells are assumed to be twice as tall as they are wide; `--cell-aspect` overrides it. When the terminal reports its size in pixels (`ws_xpixel`/`ws_ypixel`), the native CLI uses the real aspect.

### Graphics protocols

Some terminals can show real images. `--encoding kitty` uses the [Kitty graphics protocol](https://sw.kovidgoyal.net/kitty/graphics-protocol/) (raw RGBA, base64 in 4096-byte chunks), `iterm2` the [iTerm2 inline images protocol](https://iterm2.com/documentation-images.html) (a PNG), and `sixel` DEC Sixel graphics with a 256-colour palette from a deterministic median cut. They are used only when asked for — a terminal's `TERM` does not say reliably whether it supports them — and the image is sized from the cell count and the cell's pixel size.

Inside tmux the sequences have to be wrapped in tmux's passthrough, which since tmux 3.3 also needs `set -g allow-passthrough on` ([tmux FAQ](https://github.com/tmux/tmux/wiki/FAQ)). `snapshot` refuses a graphics protocol inside tmux unless `--tmux` asks for the wrapping; the cell encodings need none of this.

## Camera, shading and overlays

| Option (CLI / Python) | Default | Meaning |
| --- | --- | --- |
| `--view` / `view` | — | `iso`, `+x`, `-x`, `+y`, `-y`, `+z`, `-z`: the camera sits on that side and looks at the mesh |
| `--azimuth`, `--elevation`, `--roll` | 45, 35.26, 0 | the camera direction in degrees, the SVG writer's convention (the defaults are the isometric view) |
| `--perspective` / `projection="perspective"` | orthographic | a pinhole camera, with `--fov` (30°) vertical field of view |
| `--zoom`, `--pan-x`, `--pan-y` / `zoom`, `pan` | 1, 0, 0 | magnify the fitted view; shift it by a fraction of the frame |
| `--shading` / `shading` | `flat` | `none` (unlit), `flat` (one Lambert term per face) or `smooth` (normals from [`compute_normals`](./normals.md), split at `--split-angle`, 30°) |
| `--one-sided` / `two_sided=False` | two-sided | leave back faces at the ambient level instead of lighting both sides |
| `--ambient` / `ambient` | 0.25 | the fraction of its colour a face turned away from the light keeps |
| `light_dir` | `(0, 0, 1)` | direction toward the light in camera space; the default is a headlight |
| `--edges` / `edges` | `none` | `all` element edges, or the `feature` edges of [`feature_edges`](./feature_edges.md) (sharp at `--feature-angle`, open and non-manifold) |
| `--supersample` / `supersample` | 1 | 2 or 4 samples per pixel along each axis, averaged by a box filter |
| `--axes`, `--scale-bar` | off | the world axes (x red, y green, z blue) in the lower left; a 1-2-5 scale bar in the lower right (orthographic only) |
| `--fill`, `--line-color`, `--edge-color`, `--background` | `#c8c5bd`, `#000080`, `#000000`, `none` | colours as `#rrggbb`, `#rrggbbaa` or `none` (transparent); Python takes `(r, g, b[, a])` tuples |

The camera is fitted to the geometry actually drawn, not to the whole bounding box of the points. Shading multiplies the face colour (the fill or the mapped field) and never replaces it.

## Colouring by data

The [`--color-by` family](./formats/svg.md) of `convert` works the same way here: `--color-by NAME` (point data first, then cell data), `--component`, `--cmap`, `--vmin`, `--vmax`, `--nan-color` and `--colorbar`. Point data colours a face by the mean of its corners, cell data by the owning cell, and non-finite values take `--nan-color`. The automatic range spans every face handed to the rasterizer, hidden ones included, so an orbit never changes the colours. The range is printed under a text rendering (`--no-notes` turns that off) and returned by `render_image(..., return_info=True)`.

The colormaps are `viridis` (the default), `coolwarm`, `turbo`, `magma`, `inferno`, `plasma` and `grey`, and a reversed `_r` variant of each. They are shared with the SVG, TikZ and glTF writers. viridis, magma, inferno and plasma are matplotlib's listed maps (CC0); Turbo is Apache-2.0 (Google LLC); coolwarm and grey are sampled from matplotlib.

## Field rendering

These options (v16.34.0) turn a picture of a surface into a picture of a solution. Each is a flag of `snapshot` and a keyword argument of the same name in every language (see [surfaces](#surfaces)).

| Option (CLI / Python) | Meaning |
| --- | --- |
| `--expr TEXT` / `expr` | colour by a [`data_calc`](./data_calc.md) expression (`"mag(u) / max(p)"`), evaluated on the fly: operands are point data, then cell data, so no intermediate file is written |
| `--reduce NAME` / `reduce` | colour a six- or nine-component tensor array by `mises`, `hydrostatic` or `principal` (with `--component` 0, 1, 2 the smallest, middle and largest, the largest by default), as [`tensor_invariants`](./tensor_invariants.md) computes them |
| `--clip LOW,HIGH` / `clip=(low, high)` | bound the automatic range by percentiles of the drawn values, so one outlier does not flatten the plot; either end may be left empty (`--clip 2,`); an explicit `--vmin`/`--vmax` wins |
| `--symmetric` / `symmetric` | make the automatic range symmetric about zero, so a diverging map such as `coolwarm` puts zero at its midpoint |
| `--scale linear\|log\|symlog`, `--scale-threshold T` / `scale`, `scale_threshold` | the colormap position of a value: proportional, proportional to its logarithm (a value that is not positive is drawn in `--nan-color`, and a range that is not positive is an error), or a signed logarithm that is linear within `T` of zero |
| `--colorbar` | a gradient bar with tick marks; the tick values and the range are also in the notes under a text picture, and a count of values with no finite value on the scale says so |
| `--categorical` | treat the integer array of `--color-by` (a material id) as categories, coloured from a fixed qualitative palette, with a key in the notes |
| `--color-regions` | colour cells by the first named cell [region](./regions.md) that holds them (regions in their sorted order), with a key listing each name, colour and cell count; an empty region still appears in the key |
| `--category-edges` | draw the edges where two faces of different category meet (categories, regions or a flag diagnostic) |
| `--isolines N`, `--iso-levels A,B,...`, `--iso-color` / `isolines`, `iso_levels`, `iso_color` | contour lines of the point array `--color-by`: `N` equally spaced levels inside the range, or the explicit levels, drawn with a depth bias so they sit on their faces; cell data are refused by name |
| `--vectors NAME`, `--vector-count N`, `--vector-length L`, `--vector-color` / `vectors`, ... | arrows for a vector point array at about `N` evenly ranked drawn points; the longest is 6% of the model's diagonal, or every arrow is `L` model units long; the heads lie in the plane of the arrow and the line of sight |
| `--warp NAME`, `--warp-scale S`, `--warp-outline`, `--outline-color` / `warp`, ... | move the points by a displacement point array times `S`, optionally drawing the undeformed outline (its open and sharp edges) beside the deformed shape |
| `--diagnostic NAME`, `--quality-metric M` / `diagnostic`, `quality_metric` | one flag for the checks people run first, below |

The diagnostics reuse metrics that already exist; none is a new algorithm:

| `--diagnostic` | Draws |
| --- | --- |
| `quality` | cells coloured by the [`quality`](./mesh_quality.md) metric `--quality-metric` names (`scaled_jacobian`, `aspect_ratio`, ...) |
| `inverted`, `degenerate` | cells with a negative signed volume or area, or a near-zero volume, edge or Jacobian, in red; the notes count the faces |
| `orientation` | front faces (towards the camera) in blue and back faces in orange, so a flipped patch is visible at once |
| `free-edges` | open edges in orange, non-manifold edges in red and inconsistently wound pairs in magenta, from the [`repair`](./repair.md) detectors |
| `edge-length` | faces coloured by their mean edge length |

A diagnostic excludes `--color-by`, `--expr` and `--color-regions`; the options that make no sense together are refused by name rather than ignored. With none of the new options set a frame is byte-identical to v16.33.0's.

## Surfaces

The rasterizer is reachable from every language the library has, with the same options and the same pixels:

| Surface | Entry points | See |
| --- | --- | --- |
| native CLI, Python CLI | `snapshot` | [CLI](./cli.md#meshioplusplus-snapshot) |
| Python | `render_image`, `render_text`, `snapshot` (core only) | below |
| MCP | `render_mesh` | [MCP server](./mcp.md) |
| C | `mio_render`, `mio_frame_*`, `mio_render_text`, `mio_write_snapshot` | [C API](./c_api.md#software-rendering-v16-34-0) |
| Fortran | `mio_mesh%render`, `%render_text`, `%write_snapshot`, `mio_frame` | [Fortran](./fortran.md#software-rendering-v16-34-0) |
| Julia | `render`, `render_text`, `render_png`, `write_snapshot` | [Julia](./julia.md#software-rendering-v16-34-0) |
| R | `mio_render`, `mio_render_text`, `mio_render_png`, `mio_write_snapshot` | [R](./r.md#software-rendering-v16-34-0) |
| WebAssembly | `render`, `renderText`, `renderPng` (cell encodings only) | [WebAssembly](./wasm.md#software-rendering-v16-34-0) |
| settings pipeline | the `Snapshot` step | [pipeline](./pipeline.md) |

The interactive loop and the graphics-protocol encoders are CLI features: the library carries frames, text and PNG, never a tty. The `Snapshot` pipeline step writes a frame as a side output after any step, with the mesh passing through untouched, so a `Smooth` or `Decimate` can be shown before and after in a report.

## The interactive viewer (`tui`)

![The interactive viewer showing the bunny coloured by height, with a colour bar and a status line](/images/tui_loop.svg)

*`meshioplusplus tui bunny.stl --color-by height --colorbar --shading smooth` after a mouse drag, in a 100 x 34 terminal. This and the other figures on this page are generated by replaying a recorded input stream (`tools/gen_doc_images.py`), so they show exactly what the viewer writes.*

`tui` takes over the terminal and lets you look around a mesh with the mouse and the keyboard. It is the viewer for a machine with no display: nothing to install, no GPU, no X forwarding.

```sh
meshioplusplus tui part.vtu
meshioplusplus tui result.vtu --color-by temperature --colorbar --cmap magma --edges feature
meshioplusplus tui part.vtu --encoding braille --color-depth 256   # a plainer terminal
```

```python
import meshioplusplus as mio

mio.tui(mio.read("part.vtu"), color_by="temperature")   # returns when you quit
mio.view(mesh, backend="terminal")                      # the same viewer through view()
```

It takes the render flags of [`snapshot`](#field-rendering) (`--color-by`, `--cmap`, `--edges`, `--shading`, `--view`, …), plus `--encoding`, `--color-depth`, `--cell-aspect` and `--tmux`; the size of the picture is the size of the terminal.

### Keys and mouse

| input | does |
| --- | --- |
| drag with the left button | orbit: a drag across the whole width turns the camera 360°, down the whole height 180° |
| drag with the right button, or shift-drag | pan |
| mouse wheel, `+` `-`, Page Up / Page Down | zoom |
| arrow keys | pan |
| `1` … `7` | the views `iso`, `+x`, `-x`, `+y`, `-y`, `+z`, `-z` |
| `p` | switch between orthographic and perspective |
| `e` | edges: off, all, feature (re-prepares the scene) |
| `s` | shading: flat, smooth, none (smooth re-prepares the scene) |
| `a` `b` `c` | the axes, the scale bar and the colour bar |
| `r`, Home | reset the camera |
| `?` | show or hide the help |
| `q`, Escape, Ctrl-C, Ctrl-D | quit |

![The viewer with its help panel open](/images/tui_help.svg)

The status line shows the camera (`az`, `el`, `zoom`, the projection) and the last thing you changed; the lines above it are the notes of the frame (the colour range, ticks and keys), as in `snapshot`.

### What it does for you

- **The mesh is prepared once.** The boundary skin of a volume, the field and its range, the normals and the edge lines do not depend on the camera, so the viewer builds them when it starts ([`prepare_render`](#prepared-scenes)) and every frame after that is only the projection and the raster. At tier XL, where the skin costs seconds, this is the difference between a viewer and a slideshow. Toggling edges, or smooth shading for the first time, prepares again, and says so.
- **The model keeps its size.** The frame is fitted to the bounding sphere of what is drawn, not to the projected extents, so the object turns in place instead of growing and shrinking as it rotates.
- **Only changed cells are sent.** Each frame is compared with the last (cell by cell, as the chosen colour depth shows them) and only the differences are written, with cursor addressing, in one write. That keeps an SSH link usable.
- **The terminal is always given back.** Raw mode, the alternate screen, the hidden cursor, mouse reporting and bracketed paste are undone on a normal quit, on an error, on `SIGINT`, `SIGTERM` and `SIGHUP` (the viewer sees the signal within 50 ms and leaves cleanly, exiting with 128 plus the signal number), and at `exit`. Only `SIGKILL` cannot be undone; `reset` fixes that. A second signal while the first is being handled restores the terminal on the spot.
- **It refuses to guess.** With standard input or output not a terminal it stops with a message naming [`snapshot`](#outputs), which draws one frame to a file or a pipe, and changes nothing.
- **Resizing refits.** `SIGWINCH` (or, on Windows, a size poll) redraws at the new size; a picture is never larger than the screen, and a screen too small for a picture says so.

The final frame at rest is the deterministic one: the same camera gives the same cells whatever path led to it.

![The viewer in Braille encoding with feature edges](/images/tui_braille.svg)

### Choosing the encoding

`--encoding` is the same set as for `snapshot`. `halfblock` (the default) is the portable choice; `quadrant`, `sextant` and `braille` pack more detail per cell and need a font that has the glyphs; `ascii` works anywhere. The graphics protocols (`kitty`, `iterm2`, `sixel`) redraw the whole image on every change, so they suit a local terminal rather than a slow link. Inside tmux they need `--tmux` (and `set -g allow-passthrough on`), and the viewer refuses by name otherwise.

### Prepared scenes

The library exposes the cached half as `prepare_render(mesh, options)` and `render_scene(scene, options, fixed_fit)` in `operations/render.hpp`: `prepare_render` fixes the field, edges, colours and shading mode; `render_scene` takes the camera, the frame size, the lighting and the overlays. With the same options the frame is byte-identical to `render`'s, which a test checks across cameras, fields, vectors, warps, diagnostics and supersampling. `encode_cells` and `encode_cells_update` give the cell grid and the bytes that repaint one grid over another. They are what the viewer uses, and available to anyone writing their own loop; the loop itself, with its raw mode and signal handling, stays in the command-line layer.

### Replaying a session

For tests and documentation, `--replay FILE` plays a recorded input stream on a `--cols` by `--rows` screen instead of a terminal and prints everything the viewer wrote. `tui(mesh, replay=bytes, cols=…, rows=…)` does the same from Python and returns the output with the final camera. The figures on this page, and the viewer's own tests, are made this way: the bytes are decoded by `tools/ansi_screen.py`, a small terminal emulator that turns them into a screen to assert on, or into an SVG.

### Windows

The viewer asks the console for virtual-terminal input and output with `SetConsoleMode` (`ENABLE_VIRTUAL_TERMINAL_INPUT` and `ENABLE_VIRTUAL_TERMINAL_PROCESSING`), sets the code page to UTF-8 for the session and restores all of it on the way out; a console that refuses either flag gets a message naming [`snapshot`](#outputs). Microsoft documents [what the console emits for keys in that mode](https://learn.microsoft.com/en-us/windows/console/console-virtual-terminal-sequences#input-sequences) (the arrow keys, Home and End, Page Up and Down, Insert, Delete, the function keys, Ctrl and Alt), which is what the key parser reads. That page says nothing about mouse reports, so mouse input on Windows is not verified: Windows Terminal is expected to forward the SGR mouse mode, the legacy console host may not, and the keys and the wheel-less `+` `-` zoom work either way.

This path is built and smoke-tested in CI, but a runner has no console to drive, so it has not been exercised on a real one. Please report what you see on Windows.

## Screenshots without a viewer

[`screenshot()`](./viewer.md#screenshots) needs Polyscope and EGL or a virtual framebuffer. Without Polyscope it now falls back to this rasterizer, with a warning, and the `screenshot` verbs of both CLIs do the same. The native release binaries, which exclude Polyscope, therefore write a PNG instead of reporting a build flag, and their `view` error points at `snapshot`.

## Determinism

A frame is a pure function of the mesh and the options: the same pixels from the native CLI, Python, every parallel backend and every thread count. The rules that make it so:

- Screen positions are rounded once to fixed point with 8 subpixel bits, so coverage is decided by exact 64-bit integer edge functions. A pixel centre exactly on an edge belongs to the triangle for which that edge is a top or left edge, so two triangles sharing an edge never both claim, or both miss, a pixel.
- Depth and shading are interpolated in double precision from those integers in a fixed expression order. The whole core is compiled with `-ffp-contract=off`, so no fused multiply-add changes a rounding between compilers.
- The depth test is strict: on a tie the primitive drawn first keeps the pixel, and primitives are drawn in block-major order (triangles, then lines, then points).
- The frame is cut into 32x32-pixel tiles, each owned by one worker; there is no shared write to order.
- Supersampling is resolved by an integer box filter in a fixed order, and the text encoders break ties in a fixed order too.

The operation-determinism check (`tools/bench_ops.sh … --hash`, see [benchmarks](./benchmarks.md#determinism-check)) holds the rasterizer to this through its `render_frame_160x96` and `render_ss4_800x480` rows, which hash both the pixels and the id buffer.

## Python and MCP

`render_image(mesh, width, height, **options)` returns the `(H, W, 4)` uint8 image; `return_info=True` adds `cell_ids` (an `(H, W)` array of the input cell drawn at each pixel, global block-major, -1 for none), the mapped range and the notes. `render_text(mesh, cols, rows, encoding, color_depth, format, **options)` returns the text, sized to the terminal by default; `format` is `ansi`, `plain` or `html`. `snapshot(mesh, path, **options)` writes any of the forms above. They need the compiled core: there is no pure-Python rasterizer, by decision (see below), so without `_core` they raise `NotImplementedError`.

The [MCP server](./mcp.md)'s `render_mesh` tool lets an agent look at a mesh: without an output path it returns plain text, with a `.png` path an image. The interactive `tui` is deliberately not a tool: it hands the terminal to a person.

## Design decisions

These were settled before the rasterizer was written ([roadmap](./roadmap.md) §7.1).

**Frame budget.** Measured with the `render_frame_160x96` (a terminal-sized frame of the tier's volume) and `render_ss4_800x480` (a 4x-supersampled 800x480 image of the same volume with smooth shading, every edge and a mapped field) rows of the [benchmark harness](./benchmarks.md), on four OpenMP threads of a 16-core machine:

| tier | cells | 160x96 | 4x supersampled 800x480 |
| --- | --- | --- | --- |
| M | 162,000 | 21 ms | 50 ms |
| L | 750,000 | 103 ms | 144 ms |
| XL | 10,368,000 | 3.4 s | 3.6 s |

At XL both rows cost about the same although the second draws 400 times as many samples, so what they measure is the boundary-skin extraction, not the raster: a frame of a few tens of milliseconds is already inside the budget of an interactive loop. The consequence for the interactive loop is that the skin (and the triangle list built from it) is extracted once per mesh and cached, and that a coarser draw while the camera moves buys nothing at these sizes: the tile rasterizer needs no spatial structure, since a primitive costs only the tiles it touches. Whether a coarse draw helps on a much slower terminal link is a question about output, not rendering, and belongs to the loop.

**No Python twin.** Every other renderer in meshio++ (SVG, TikZ, glTF) has a byte-pinned NumPy twin. This one does not: a NumPy edge-function rasterizer would be slow and would double the surface that has to stay byte-identical. Python calls the core or raises by name, the precedent of [`remesh_volume`](./remesh_volume.md) and the MMG operations; `MESHIOPLUSPLUS_STRICT_CORE` and the core-fallback guards are unaffected because there is no fallback to take.

**Nothing beyond the standard library and the operating system.** No curses, image or terminal library. The PNG encoder, the colour quantizers, the Sixel median cut and the asciicast writer are all part of the core. The terminal itself is queried only by the CLI layer, through `ioctl(TIOCGWINSZ)` on POSIX and the console API on Windows. On Windows the CLI turns on `ENABLE_VIRTUAL_TERMINAL_PROCESSING` and the UTF-8 output code page; [virtual-terminal processing](https://learn.microsoft.com/en-us/windows/console/console-virtual-terminal-sequences) needs Windows 10 version 1511 (build 10586) or later, and the interactive loop also needs virtual-terminal *input*, for which Microsoft's documentation gives no version (the viewer asks and, if the console refuses, says so). On an older console `snapshot -` falls back to glyphs without colour.

**Where the code lives.** The public header is `operations/render.hpp` (`RenderOptions`, `Frame`, `TextOptions`, `render`, `encode_text`, `encode_png`, `write_snapshot`, …); the rasterizer, PNG encoder and text encoders are core-private (`src/cpp/src/detail/raster.*`, `png_write.*`, `operations/render_text.cpp`). Terminal queries live in the CLI (`src/cpp/cli/terminal.*`), so the library never owns a tty. The header is an additive ABI change, recorded in the [ABI reviews](./abi_reviews.md).

**Determinism.** As above: 8 subpixel bits, integer edge functions with a top-left fill rule, a strict depth test, tile ownership and an integer box filter. A separate `-ffp-contract=off` for the rasterizer turned out to be unnecessary: the core is already built with it.

## Limits

- No transparency, shadows or reflections; one directional light.
- The scale bar is drawn for the orthographic camera only.
- Screen coordinates are clamped to ±262,144 pixels, so an extreme zoom distorts geometry far outside the frame (never inside it).
- Streamlines are not drawn: they need a sampler and a step integrator. Vector arrows sit at the drawn points, not along a flow.
- Isolines are the piecewise-linear contours of the interpolant on each triangle (a quad is split on its shorter diagonal); they are not smoothed.
- The id buffer of a supersampled frame reports the cell at the centre sample of each pixel.
- The interactive viewer has no clip planes, probing, command line, sessions, comparison, time stepping or `--follow` yet; the [roadmap](./roadmap.md) lists them.
- The Kitty, iTerm2 and Sixel encoders are written to their specifications but are not tested against every terminal; the cell encodings are the portable choice.
