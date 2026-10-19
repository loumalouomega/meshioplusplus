# Third-party components in the meshio++ binaries

meshio++ itself is MIT ([LICENSE](LICENSE)). The binaries the project publishes — the PyPI wheels (the `meshioplusplus._core` extension) and the release command-line executables — are statically linked with the libraries below. A build from source links whichever of the optional native libraries it finds (see [installation](doc/installation.md)); the WebAssembly package lists its own in [`src/wasm/THIRD_PARTY_NOTICES.md`](src/wasm/THIRD_PARTY_NOTICES.md).

| Component | Version | Licence | Wheels | CLI binaries |
| --- | --- | --- | --- | --- |
| [Open CASCADE Technology](https://github.com/Open-Cascade-SAS/OCCT) | 7.6.3 | LGPL-2.1-only WITH OCCT-exception-1.0 | yes | yes |
| [zlib](https://zlib.net/) | 1.3.1 | Zlib | yes | yes |
| [gidpost](https://www.gidsimulation.com/) (vendored) | 2.14 | BSD-2-Clause | yes | yes |

**OpenCASCADE is LGPL-2.1 with the OCCT exception**, so the published binaries are not plain MIT: their licence is `MIT AND (LGPL-2.1-only WITH OCCT-exception-1.0)`. The meshio++ sources, and a build from source that does not link OpenCASCADE (`-DMESHIOPLUSPLUS_WITH_OCCT=OFF`, or no OpenCASCADE found), are MIT.

## OpenCASCADE and the LGPL

OpenCASCADE's licence is the GNU Lesser General Public License 2.1 (`LICENSE_LGPL_21.txt`) with the additional exception in `OCCT_LGPL_EXCEPTION.txt`, which permits using its headers' inline code and templates without the resulting object code becoming a "work based on the library". Both texts ship beside this notice in every wheel (its `.dist-info/licenses`) and every CLI archive.

The binaries link OpenCASCADE statically, as the CAD kernel of `fit_cad` and of STEP and IGES reading. To let a recipient exercise the LGPL's right to relink against a modified OpenCASCADE, every binary is reproducible from source:

1. The OpenCASCADE source is the upstream tag `V7_6_3`, commit `b079fb9877ef64d4a8158a60fa157f59b096debb`, unmodified (no patches are applied).
2. `build/build-occt.sh PREFIX` fetches that commit, verifies it, and builds the static libraries (`build/occt-recipe.sh` holds the pin and the configuration). Configuring meshio++ with `OpenCASCADE_ROOT=PREFIX` (or `-DOpenCASCADE_DIR=PREFIX/lib/cmake/opencascade`) links them: `pip install .` for the wheel, `build/configure.sh --cli --build` for the executable. The workflows that build the published binaries, `.github/workflows/wheels.yml` and `.github/workflows/cli.yml`, run exactly these steps. All of it is in the meshio++ repository at the tag of the binary's version.
3. Replacing the OpenCASCADE checkout with a modified one (or pointing `OpenCASCADE_ROOT` at any compatible OpenCASCADE 7.6 or later, static or shared) and rebuilding produces a binary linked against it.

## Licence texts

OpenCASCADE's are named above. zlib's is in the header of its source distribution at the version above (the pin is in `CMakeLists.txt`, `MESHIOPLUSPLUS_ZLIB_STATIC`); the vendored gidpost's is `src/cpp/third_party/gidpost/LICENSE` in the meshio++ repository.
