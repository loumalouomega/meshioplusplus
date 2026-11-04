# The pinned upstream versions of the libraries the wheels and the WebAssembly
# artifacts link, sourced by build/build-wasm-deps.sh (Emscripten) and
# build/build-permissive-libs.sh (the native static prefix of the PyPI wheels)
# so the two cannot drift. SHA256s are of the release tarballs; the commits are
# what a tag peels to, for a clone where an archive download is refused.
# HDF5 >= 1.14.0 is required, not merely preferred: it removed
# H5detect/H5make_libsettings.
HDF5_VERSION="1.14.6"
NETCDF_VERSION="4.9.3"
CGNS_VERSION="4.5.2"
HDF5_SHA256="e4defbac30f50d64e1556374aa49e574417c9e72c6b1de7a4ff88c4b1bea6e9b"
NETCDF_SHA256="990f46d49525d6ab5dc4249f8684c6deeaf54de6fec63a187e9fb382cc0ffdff"
CGNS_SHA256="95075e1fd0b51d97b1b96b73ebe03b1a551fbcc9cd2b2b6f487ccccedcff5964"
NETCDF_COMMIT="3a6b1debf1557f07b606ce3653e44f0d711203be"
CGNS_COMMIT="3a972d92559787a49c5f9a54974ea1946ecb9d21"
ZSTD_VERSION="1.5.7"
ZSTD_SHA256="eb33e51f49a15e023950cd7825ca74a4a2b43db8354825ac24fc1b7ee09e6fa3"
LZ4_VERSION="1.10.0"
LZ4_SHA256="537512904744b35e232912055ccf8ec66d768639ff3abe5788d90d792ec5f48b"
KAHIP_VERSION="3.25"
KAHIP_SHA256="3abad20158887b585b4d4792c35fa9023b7fd634b2bc4494334a3951b69f4842"
KAHIP_COMMIT="2e86d3da2065bcf6ddbf6fb84d2b52f33db1061f"
NETGEN_VERSION="6.2.2401"
NETGEN_COMMIT="2ff886457fe899685d9d826c8a29f58a31889dda"
# MMG (the full WebAssembly artifact's adapt, discretize_levelset and
# mesh_domain): a clone of its tag, refused unless it is MMG_COMMIT, as Netgen
# is -- the same 5.8.0 build/build-mmg.sh builds natively.
MMG_VERSION="5.8.0"
MMG_COMMIT="4d8232c8aebfed877935d75d4d4a67e850962422"
# QuadriFlow (MIT), the quad remesher remesh_quads runs as a separate program
# (roadmap 7.7.3). Never linked or shipped: CI builds the executable from this
# commit to test the backend. It has no release tags.
QUADRIFLOW_COMMIT="810b7a0967c35b0dc85b4464e3835e26a756c967"
# CutCells (MIT), the backend of cut_level_set and cut_quadrature (doc/cut.md):
# a clone of this commit (the release-0.5.0 merge; CutCells has no release
# tags), built static by build/build-cutcells.sh for the wheels and CI.
CUTCELLS_VERSION="0.5.0"
CUTCELLS_COMMIT="3ccd9b03a1b2a607f137fc554057ccd002e4bc9b"
