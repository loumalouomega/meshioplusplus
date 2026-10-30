# Licence of the published binaries

meshio++ itself is MIT ([LICENSE](LICENSE)). The base PyPI wheel, `meshioplusplus` (the `meshioplusplus._core` extension and the shared core beside it), and the Blender extension zip that carries it, contain no copyleft code. They are `MIT AND MPL-2.0`: MPL-2.0 is Eigen, a header-only library compiled into the Linux and macOS wheels, whose terms apply to its own files only. The other libraries the base wheel links (HDF5, netCDF, cgnslib, zstd, LZ4, bzip2, KaHIP, zlib) are permissive and listed in `THIRD_PARTY_NOTICES.md`. Other binaries the project publishes carry copyleft code: the add-on wheels `meshioplusplus-occt` (OpenCASCADE Technology 7.6.3 as shared libraries beside its extension, `MIT AND (LGPL-2.1-only WITH OCCT-exception-1.0)`), `meshioplusplus-netgen` (Netgen's nglib, `MIT AND LGPL-2.1-only`) and `meshioplusplus-mmg` (libmmg, built without ELAS, `MIT AND LGPL-3.0-or-later`), the release command-line executables, and the Conan and vcpkg packages built with their default OpenCASCADE option (these link OpenCASCADE statically). Their licence is

    MIT AND (LGPL-2.1-only WITH OCCT-exception-1.0)

OpenCASCADE is licensed under the GNU Lesser General Public License 2.1 with the additional exception in `OCCT_LGPL_EXCEPTION.txt`. Its licence texts, `LICENSE_LGPL_21.txt` and `OCCT_LGPL_EXCEPTION.txt`, ship beside this file in each of those artifacts. `pip install meshioplusplus[occt]` (and `[all]`) installs the add-on, so an environment with it holds the combined licence.

`THIRD_PARTY_NOTICES.md` lists every component each binary contains, and explains how to rebuild a binary from source against a modified OpenCASCADE.

A build from source without OpenCASCADE (`-DMESHIOPLUSPLUS_WITH_OCCT=OFF`, or the Conan option `with_occt=False`) is MIT.
