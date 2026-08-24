# Compression codecs (zlib, zstd, lz4)

VTK XML formats (`.vtu`, `.vtp`) store binary `DataArray` payloads in fixed 32 KiB blocks, each independently compressed. meshio++ supports three codecs for those blocks.

**zlib is the default and always available.** zstd and lz4 are optional, off by default, and change nothing about existing files.

## Portability — read this before choosing

| Codec | `compressor=` attribute | ParaView / VTK can read it? |
| --- | --- | --- |
| `zlib` | `vtkZLibDataCompressor` | **yes** — the default, universally supported |
| `lz4` | `vtkLZ4DataCompressor` | **yes** — a real VTK compressor |
| `zstd` | `vtkZSTDDataCompressor` | **no** — a meshio++ extension |
| `lzma` | `vtkLZMADataCompressor` | read by the Python path only |

VTK ships `vtkZLibDataCompressor`, `vtkLZ4DataCompressor` and `vtkLZMADataCompressor` — there is **no ZSTD compressor in VTK or ParaView**. (Checked against VTK 9.6: `vtkLZ4DataCompressor` exists, `vtkZSTDDataCompressor` does not.) meshio++ writes `vtkZSTDDataCompressor` by analogy so that a future VTK addition would align, and so ParaView fails with a clean "unknown compressor" rather than silently misreading. If a file has to open in ParaView, use `zlib` or `lz4`.

lz4 uses LZ4's raw **block** format (not the frame format), matching `vtkLZ4DataCompressor`. Its blocks carry no size field, so the decompressed size is taken from VTU's own block header — which is exactly what makes the raw format usable here.

## Building with them

```bash
cmake -DMESHIOPLUSPLUS_WITH_ZSTD=ON -DMESHIOPLUSPLUS_WITH_LZ4=ON ...
```

Both default **OFF**, unlike zlib/HDF5/netCDF: zlib stays the write default, so a build without them reads and writes exactly what it always did.

```python
from meshioplusplus import _core
_core.__has_zlib__, _core.__has_zstd__, _core.__has_lz4__
```

Package managers: Conan `with_zstd` / `with_lz4`, vcpkg features `zstd` / `lz4` — both off by default there too, so existing package IDs are unaffected.

**WASM**: compiled out. There is no Emscripten port for either, and zlib (`-sUSE_ZLIB=1`) is unchanged.

## Using them

```python
meshioplusplus.write("mesh.vtu", mesh)                      # zlib (default, unchanged)
meshioplusplus.write("mesh.vtu", mesh, compression="lz4")
meshioplusplus.write("mesh.vtu", mesh, compression="zstd")
```

```bash
meshioplusplus compress --codec lz4 mesh.vtu
```

`--codec` is **rejected** for formats without a block codec rather than silently ignored — believing you got zstd and quietly getting gzip would be the worst outcome:

| `--codec` | vtu / vtp | cgns / h5m / xdmf | ansys / gmsh / ply / stl / vtk |
| --- | --- | --- | --- |
| *unset* | today's behaviour | today's behaviour | today's behaviour |
| `zlib` | zlib | accepted (gzip already is zlib) | **error** |
| `lz4` / `zstd` | new path | **error** | **error** |

With no `--codec`, every format behaves exactly as before.

## Reading a file whose codec you lack

The codec is detected from the file's `compressor=` attribute. If this build cannot handle it, the error names the option to enable:

```
VTK XML lz4 decompression requires a build with -DMESHIOPLUSPLUS_WITH_LZ4=ON
```

The pure-Python reference path also supports both codecs, via optional dependencies:

```bash
pip install "meshioplusplus[codecs]"     # lz4 + zstandard
```

A file needing a codec present in **neither** the C++ core nor Python is genuinely unreadable. That is a new failure class introduced by this feature, and it fails by name rather than with an `ImportError` from inside the parser.

### Interop is tested, not assumed

`tests/python/test_codecs.py` runs against a real VTK when one is installed (`pip install vtk`) and asserts all three directions:

- VTK reads a `.vtu` meshio++ wrote with `lz4`;
- meshio++ reads a `.vtu` VTK wrote with `SetCompressorTypeToLZ4()`;
- VTK **cleanly refuses** a `zstd` file (`Error creating vtkZSTDDataCompressor`) rather than misreading it — which is the trade being made by writing a non-VTK compressor name.

## Compression is framing-only

All codecs share VTU's block framing verbatim, and every codec's decoded payload is required by test to be byte-identical to the uncompressed variant's. Compression changes how the bytes are packaged, never which bytes they are — files written by the C++ core and by the pure-Python reference read interchangeably in both directions.
