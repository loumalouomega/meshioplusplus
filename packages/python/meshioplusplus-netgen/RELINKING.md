# Relinking and source for Netgen's nglib

This wheel (`meshioplusplus-netgen`) contains Netgen's nglib and ngcore, licensed LGPL-2.1 (`COPYING.LESSER`). The meshio++ extension is MIT (`LICENSE`).

- nglib is linked **dynamically**. Its shared libraries ship in this wheel beside the extension, and the extension finds them through its RPATH, so you can replace them with your own build of a compatible Netgen without rebuilding the extension.
- Version: Netgen 6.2.2401, commit `2ff886457fe899685d9d826c8a29f58a31889dda` (`NETGEN_VERSION` and `NETGEN_COMMIT` in `build/deps-pins.sh`).
- Source: https://github.com/NGSolve/netgen at that commit, unmodified (the patch in `build/patches/` is for the WebAssembly build only).
- Build recipe: `build/build-netgen.sh` in the meshio++ repository (https://github.com/loumalouomega/meshioplusplus), at the tag matching this wheel's version.
- To rebuild the extension: install the matching `meshioplusplus` wheel and build `packages/python/meshioplusplus-netgen` from the repository with `NETGEN_ROOT` pointing at your Netgen install.
