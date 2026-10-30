"""Smoke test of the base wheel alone (wheels.yml's CIBW_TEST_COMMAND;
standalone, no pytest): it carries no OpenCASCADE, so STEP reading raises an
error naming the add-on, ``fit_cad`` falls back to the native backend, the
shared core and its CMake package are in the wheel, and (Linux and macOS) the
native HDF5, netCDF, cgnslib and codec paths are active without their symbols
leaking, beside h5py and netCDF4 in the same process."""

import os
import re
import shutil
import subprocess
import sys
import tempfile

import numpy as np

import meshioplusplus as mp
from meshioplusplus import _core
from meshioplusplus._fit_cad import has_occt

assert (
    not has_occt() and not _core.__has_occt__
), "the base wheel must carry no OpenCASCADE"

path = os.path.join(tempfile.mkdtemp(), "cube.step")
points = np.array(
    [[x, y, z] for z in (0.0, 1.0) for y in (0.0, 1.0) for x in (0.0, 1.0)]
)
quads = np.array(
    [[0, 2, 3, 1], [4, 5, 7, 6], [0, 1, 5, 4], [1, 3, 7, 5], [3, 2, 6, 7], [2, 0, 4, 6]]
)
mp.write(path, mp.Mesh(points, [("quad", quads)]))
try:
    mp.read(path)
except NotImplementedError as exc:
    assert "meshioplusplus-occt" in str(exc), exc
else:
    raise AssertionError("STEP reading must need the add-on")

cmake_dir = mp.get_cmake_dir()
assert os.path.isfile(os.path.join(cmake_dir, "meshioplusplusConfig.cmake")), cmake_dir

# --- the permissive libraries (roadmap 4.7.2): Linux and macOS wheels ------------
# Windows keeps HDF5, netCDF and the codecs on h5py, netCDF4 and the stdlib.
if sys.platform != "win32":
    for flag in (
        "__has_hdf5__",
        "__has_netcdf__",
        "__has_cgnslib__",
        "__has_zstd__",
        "__has_lz4__",
        "__has_bzip2__",
        "__has_kahip__",
    ):
        assert getattr(_core, flag, False), f"the wheel was built without {flag}"

    # None of those libraries' symbols may leave the shared core: h5py and netCDF4
    # bring their own HDF5 into this process, and an exported H5* would interpose
    # on it. (macOS binds each symbol to its image, so only ELF can leak.)
    if sys.platform.startswith("linux") and shutil.which("nm"):
        lib_dir = os.path.join(os.path.dirname(mp.__file__), "lib")
        core_lib = next(
            os.path.join(lib_dir, f)
            for f in sorted(os.listdir(lib_dir))
            if f.startswith("libmeshioplusplus_core") and f.endswith(".so")
        )
        exported = subprocess.run(
            ["nm", "-D", "--defined-only", core_lib],
            capture_output=True,
            text=True,
            check=True,
        ).stdout
        leaked = sorted(
            {
                line.split()[-1]
                for line in exported.splitlines()
                if re.match(
                    r"(H5|nc_|NC_|cg_|cgi_|ZSTD_|LZ4_|BZ2_|kaffpa)", line.split()[-1]
                )
            }
        )
        assert not leaked, f"the core exports its dependencies' symbols: {leaked[:10]}"

    # h5py (its own HDF5) and the native core (the static one) in one process, on
    # the same file: each reads what the other wrote. STRICT_CORE makes a fallback
    # to the h5py implementation an error, so this proves the native path.
    tri = mp.Mesh(
        np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]]),
        [("triangle", np.array([[0, 1, 2]]))],
    )
    try:
        import h5py
    except ImportError:  # pragma: no cover - the wheel tests install it
        h5py = None
    if h5py is not None:
        os.environ["MESHIOPLUSPLUS_STRICT_CORE"] = "1"
        hmf = os.path.join(tempfile.mkdtemp(), "interop.hmf")
        mp.write(hmf, tri)  # the core's HDF5
        with h5py.File(hmf, "r") as f:  # h5py's HDF5, same process, same file
            assert len(list(f.keys())) > 0
        back = mp.read(hmf)  # the core reads it again after h5py opened it
        assert back.points.shape == (3, 3)
        with h5py.File(hmf, "r+") as f:  # h5py writes, the core reads
            f.attrs["written_by_h5py"] = 1
        assert mp.read(hmf).points.shape == (3, 3)
        del os.environ["MESHIOPLUSPLUS_STRICT_CORE"]
    try:
        import netCDF4
    except ImportError:  # pragma: no cover
        netCDF4 = None
    if netCDF4 is not None:
        exo = os.path.join(tempfile.mkdtemp(), "interop.exo")
        mp.write(exo, tri)  # the core's netCDF over its HDF5
        with netCDF4.Dataset(exo) as ds:  # netCDF4's own copies of both
            assert "num_nodes" in ds.dimensions
        assert mp.read(exo).points.shape == (3, 3)

print("base wheel smoke OK:", mp.__version__)
