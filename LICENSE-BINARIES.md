# Licence of the published binaries

meshio++ itself is MIT ([LICENSE](LICENSE)). The base PyPI wheel, `meshioplusplus` (the `meshioplusplus._core` extension and the shared core beside it), and the Blender extension zip that carries it, are MIT: they contain no copyleft code. Other binaries the project publishes are not MIT alone: the `meshioplusplus-occt` add-on wheel (it ships OpenCASCADE Technology 7.6.3 as shared libraries beside its extension), the release command-line executables, and the Conan and vcpkg packages built with their default OpenCASCADE option (these link OpenCASCADE statically). Their licence is

    MIT AND (LGPL-2.1-only WITH OCCT-exception-1.0)

OpenCASCADE is licensed under the GNU Lesser General Public License 2.1 with the additional exception in `OCCT_LGPL_EXCEPTION.txt`. Its licence texts, `LICENSE_LGPL_21.txt` and `OCCT_LGPL_EXCEPTION.txt`, ship beside this file in each of those artifacts. `pip install meshioplusplus[occt]` (and `[all]`) installs the add-on, so an environment with it holds the combined licence.

`THIRD_PARTY_NOTICES.md` lists every component each binary contains, and explains how to rebuild a binary from source against a modified OpenCASCADE.

A build from source without OpenCASCADE (`-DMESHIOPLUSPLUS_WITH_OCCT=OFF`, or the Conan option `with_occt=False`) is MIT.
