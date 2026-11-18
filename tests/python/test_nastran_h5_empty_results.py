"""A Nastran HDF5 result table with no row for a block must still be written.

A mesh whose per-cell result exists for one block only (NaN for the other) leaves
the other block's table empty; the Python writer reshaped the empty columns with
``reshape(0, -1)``, which NumPy cannot infer (found by the notebook job, roadmap
3.11.4).
"""

import numpy as np
import pytest

import meshioplusplus as mp

h5py = pytest.importorskip("h5py")


def _mesh():
    ply = np.array([[120.0, 45.0, -30.0, -95.0]])
    return mp.Mesh(
        np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [2, 0, 0]], dtype=float),
        [("quad", np.array([[0, 1, 2, 3]])), ("line", np.array([[1, 4]]))],
        cell_data={
            "nastran:eid": [np.array([1]), np.array([2])],
            "nastran:pid": [np.array([5]), np.array([6])],
            "nastran:card": [np.array([48]), np.array([24])],
            "nastran:x": [np.zeros((1, 3)), np.array([[0.0, 1.0, 0.0]])],
            "STRESS:X1@ply": [ply, np.full((1, 4), np.nan)],
            "STRESS:X1": [np.array([ply[0, 0]]), np.array([np.nan])],
            "FORCE:AXIAL@station": [np.full((1, 2), np.nan), np.array([[10.0, -4.0]])],
        },
    )


def test_a_result_present_on_one_block_only_is_written(tmp_path):
    path = tmp_path / "mixed.h5"
    mp.write(path, _mesh(), file_format="nastran_h5")
    back = mp.read(path)
    assert sorted(b.type for b in back.cells) == ["line", "quad"]
    assert "STRESS:X1@ply" in back.cell_data
