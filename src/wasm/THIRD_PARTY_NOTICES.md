# Third-party components in @meshioplusplus/wasm

meshio++ itself is MIT ([LICENSE](https://github.com/loumalouomega/meshioplusplus/blob/master/LICENSE)). The WebAssembly artifacts in `dist/` are statically linked with the libraries below. The package ships two pairs of artifacts, and which libraries an application carries depends on which pair it loads:

- **Lean** (`meshioplusplus_wasm`, `meshioplusplus_wasm_mt`; `loadMeshioPlusPlus()`, the default): every component here is under a permissive licence except Eigen (MPL-2.0), which the lean pair compiles in for the MED transpose; an application that loads only the lean pair takes on the MPL-2.0 file-level obligations for Eigen's own source.
- **Full** (`meshioplusplus_wasm_full`, `meshioplusplus_wasm_full_mt`; `loadMeshioPlusPlus({}, {full: true})`): everything in the lean pair, plus OpenCASCADE, zstd, lz4, bzip2, KaHIP and Netgen. **OpenCASCADE is LGPL-2.1 with the OCCT exception and Netgen is LGPL-2.1**, so the full pair is not plain MIT. The package's licence expression, `MIT AND (LGPL-2.1-only WITH OCCT-exception-1.0) AND LGPL-2.1-only`, says so.

| Component | Version | Licence | In the lean pair | In the full pair |
| --- | --- | --- | --- | --- |
| [zlib](https://zlib.net/) (Emscripten's port) | the SDK's | Zlib | yes | yes |
| [HDF5](https://www.hdfgroup.org/solutions/hdf5/) | 1.14.6 | BSD-3-Clause-style (HDF5) | yes | yes |
| [netCDF-C](https://github.com/Unidata/netcdf-c) | 4.9.3 | BSD-3-Clause | yes | yes |
| [CGNS](https://github.com/CGNS/CGNS) | 4.5.2 | Zlib | yes | yes |
| [gidpost](https://www.gidsimulation.com/) (vendored) | 2.14 | BSD-2-Clause-Views | yes | yes |
| [Eigen](https://eigen.tuxfamily.org/) (submodule) | the submodule's pinned commit | MPL-2.0 | yes | yes |
| [Open CASCADE Technology](https://github.com/Open-Cascade-SAS/OCCT) | 7.6.3 | LGPL-2.1-only WITH OCCT-exception-1.0 | no | yes |
| [Zstandard](https://github.com/facebook/zstd) | 1.5.7 | BSD-3-Clause (dual-licensed with GPL-2.0; used under BSD) | no | yes |
| [LZ4](https://github.com/lz4/lz4) (the library) | 1.10.0 | BSD-2-Clause | no | yes |
| [bzip2](https://sourceware.org/bzip2/) (Emscripten's port) | 1.0.6 | bzip2-1.0.6 | no | yes |
| [KaHIP](https://github.com/KaHIP/KaHIP) (serial kaffpa only) | 3.25 | MIT | no | yes |
| [Netgen](https://github.com/NGSolve/netgen) (nglib and ngcore only) | 6.2.2401 | LGPL-2.1-only | no | yes |

## OpenCASCADE and the LGPL

OpenCASCADE's licence is the GNU Lesser General Public License 2.1 with the additional exception in its `OCCT_LGPL_EXCEPTION.txt`, which permits using its headers' inline code and templates without the resulting object code becoming a "work based on the library". The full artifact links OpenCASCADE statically into the `.wasm` binary. To let a recipient exercise the LGPL's right to relink against a modified OpenCASCADE, the full artifact is reproducible from source:

1. The OpenCASCADE source is the upstream tag `V7_6_3`, commit `b079fb9877ef64d4a8158a60fa157f59b096debb`, unmodified (no patches are applied; it is compiled with `OCC_CONVERT_SIGNALS` undefined, since WebAssembly raises no signals to convert).
2. `build/build-wasm-deps.sh --full` fetches that commit, verifies it, and builds the wasm32 static libraries; `build/configure-wasm.sh --full-only --build` links them into `meshioplusplus_wasm_full[_mt]`. Both scripts are in the meshio++ repository at the tag of this package's version, as is every other source the artifact is built from.
3. Replacing the OpenCASCADE checkout with a modified one and re-running the two scripts produces a full artifact linked against it.

## Netgen and the LGPL

Netgen's licence is the GNU Lesser General Public License 2.1. The full artifact links its `nglib` and `ngcore` statically into the `.wasm` binary, as the backend of `tetrahedralize`; the relinking right is served the same way as OpenCASCADE's:

1. The Netgen source is the upstream tag `v6.2.2401`, commit `2ff886457fe899685d9d826c8a29f58a31889dda`, with one patch, `build/patches/netgen-6.2.2401-emscripten.patch`, which drops the JavaScript exception catching and web-only link settings its CMake imposes under Emscripten (the artifact uses WebAssembly exceptions). It is built without its GUI, Python, OpenCASCADE, MPI and CGNS interfaces.
2. `build/build-wasm-deps.sh --full` fetches that commit, verifies it, applies the patch and builds the wasm32 static libraries; `build/configure-wasm.sh --full-only --build` links them.
3. Replacing the Netgen checkout with a modified one and re-running the two scripts produces a full artifact linked against it.

An application that loads only the lean pair never fetches the full artifact's bytes, and need not ship it.

## Licence texts

Each component's licence is in its own source distribution at the version above (the pins are in `build/build-wasm-deps.sh`, Netgen's `LICENSE` at its tag; zlib and bzip2 are the ports of the pinned Emscripten SDK); the vendored gidpost's is `src/cpp/third_party/gidpost/LICENSE` in the meshio++ repository.
