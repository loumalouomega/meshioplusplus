# Licence of the published binaries

meshio++ itself is MIT ([LICENSE](LICENSE)). Some binaries the project publishes are not MIT alone: the PyPI wheels (the `meshioplusplus._core` extension), the release command-line executables, the Blender extension zip (which carries the wheel), and the Conan and vcpkg packages built with their default OpenCASCADE option. All of them link OpenCASCADE Technology 7.6.3 statically, so their licence is

    MIT AND (LGPL-2.1-only WITH OCCT-exception-1.0)

OpenCASCADE is licensed under the GNU Lesser General Public License 2.1 with the additional exception in `OCCT_LGPL_EXCEPTION.txt`. Its licence texts, `LICENSE_LGPL_21.txt` and `OCCT_LGPL_EXCEPTION.txt`, ship beside this file in each artifact.

`THIRD_PARTY_NOTICES.md` lists every component each binary contains, and explains how to rebuild a binary from source against a modified OpenCASCADE.

A build from source without OpenCASCADE (`-DMESHIOPLUSPLUS_WITH_OCCT=OFF`, or the Conan option `with_occt=False`) is MIT.
