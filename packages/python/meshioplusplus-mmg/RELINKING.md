# Relinking and source for libmmg

This wheel (`meshioplusplus-mmg`) contains MMG (libmmg, libmmg2d, libmmg3d, libmmgs), licensed LGPL-3.0-or-later (`COPYING.LESSER`, with `COPYING`, the GPL-3.0 it incorporates). The meshio++ extension is MIT (`LICENSE`).

- MMG is linked **dynamically**. Its shared libraries ship in this wheel beside the extension, and the extension finds them through its RPATH, so you can replace them with your own build of a compatible MMG without rebuilding the extension.
- Version: MMG 5.8.0, commit `4d8232c8aebfed877935d75d4d4a67e850962422` (`MMG_VERSION` and `MMG_COMMIT` in `build/deps-pins.sh`).
- Source: https://github.com/MmgTools/mmg at that commit, unmodified.
- Build recipe: `build/build-mmg.sh` in the meshio++ repository (https://github.com/loumalouomega/meshioplusplus), at the tag matching this wheel's version. It builds MMG as shared libraries **without ELAS** (GPL-3.0), so `move` is unavailable.
- To rebuild the extension: install the matching `meshioplusplus` wheel and build `packages/python/meshioplusplus-mmg` from the repository with `mmg_DIR` pointing at your MMG install.
