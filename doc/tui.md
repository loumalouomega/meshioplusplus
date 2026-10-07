# Terminal rendering

meshio++ can draw a mesh where no window, browser or GPU is available: an SSH session on a cluster login node, a container, a CI log, a terminal with nothing else installed. A software rasterizer turns the mesh into pixels, and the pixels become Unicode block characters with terminal colours, an HTML page, an asciinema recording or a PNG.

```sh
meshioplusplus snapshot part.vtu -                          # draw it in this terminal
meshioplusplus snapshot part.vtu part.png --edges feature   # a PNG, no display needed
meshioplusplus snapshot result.vtu - --color-by temperature --colorbar --cmap turbo
```

```python
import meshioplusplus as mio

mesh = mio.read("part.vtu")
print(mio.render_text(mesh, color_by="temperature"))      # ANSI text for this terminal
image = mio.render_image(mesh, 800, 600, supersample=2)   # an (H, W, 4) uint8 array
mio.snapshot(mesh, "part.png", edges="all")               # .png, .txt, .ansi, .html or .cast
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

The [MCP server](./mcp.md)'s `render_mesh` tool lets an agent look at a mesh: without an output path it returns plain text, with a `.png` path an image.

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

**Nothing beyond the standard library and the operating system.** No curses, image or terminal library. The PNG encoder, the colour quantizers, the Sixel median cut and the asciicast writer are all part of the core. The terminal itself is queried only by the CLI layer, through `ioctl(TIOCGWINSZ)` on POSIX and the console API on Windows. On Windows the CLI turns on `ENABLE_VIRTUAL_TERMINAL_PROCESSING` and the UTF-8 output code page; [virtual-terminal processing](https://learn.microsoft.com/en-us/windows/console/console-virtual-terminal-sequences) needs Windows 10 version 1511 (build 10586) or later, and the interactive loop, which also needs virtual-terminal *input*, Windows 10 1809. On an older console `snapshot -` falls back to glyphs without colour.

**Where the code lives.** The public header is `operations/render.hpp` (`RenderOptions`, `Frame`, `TextOptions`, `render`, `encode_text`, `encode_png`, `write_snapshot`, …); the rasterizer, PNG encoder and text encoders are core-private (`src/cpp/src/detail/raster.*`, `png_write.*`, `operations/render_text.cpp`). Terminal queries live in the CLI (`src/cpp/cli/terminal.*`), so the library never owns a tty. The header is an additive ABI change, recorded in the [ABI reviews](./abi_reviews.md).

**Determinism.** As above: 8 subpixel bits, integer edge functions with a top-left fill rule, a strict depth test, tile ownership and an integer box filter. A separate `-ffp-contract=off` for the rasterizer turned out to be unnecessary: the core is already built with it.

## Limits

- No transparency, shadows or reflections; one directional light.
- The scale bar is drawn for the orthographic camera only.
- Screen coordinates are clamped to ±262,144 pixels, so an extreme zoom distorts geometry far outside the frame (never inside it).
- The id buffer of a supersampled frame reports the cell at the centre sample of each pixel.
- The Kitty, iTerm2 and Sixel encoders are written to their specifications but are not tested against every terminal; the cell encodings are the portable choice.
