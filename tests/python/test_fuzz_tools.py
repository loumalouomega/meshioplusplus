"""Tests for corpus packaging and the opt-in library campaign."""

import importlib.util
import os
import pathlib
import subprocess
import zipfile

import pytest

REPO = pathlib.Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location(
    "package_seeds", REPO / "tools/fuzz/oss-fuzz/package_seeds.py"
)
seeds = importlib.util.module_from_spec(spec)
spec.loader.exec_module(seeds)


def test_packages_only_matching_small_seed_files(tmp_path):
    source, output = tmp_path / "seeds", tmp_path / "out"
    (source / "vtk").mkdir(parents=True)
    output.mkdir()
    (source / "vtk/one").write_bytes(b"valid input")
    (source / "vtk/empty").touch()
    (source / "vtk/large").write_bytes(b"x" * 262145)
    (output / "meshioplusplus_fuzz_read_vtk").touch()
    (output / "meshioplusplus_fuzz_read_vtk.dict").touch()
    seeds.package(source, output)
    with zipfile.ZipFile(
        output / "meshioplusplus_fuzz_read_vtk_seed_corpus.zip"
    ) as archive:
        assert archive.namelist() == ["one"]
        assert archive.read("one") == b"valid input"


def test_lfs_pointer_is_not_packaged(tmp_path):
    source, output = tmp_path / "seeds", tmp_path / "out"
    (source / "vtk").mkdir(parents=True)
    output.mkdir()
    (source / "vtk/pointer").write_bytes(
        b"version https://git-lfs.github.com/spec/v1\n"
    )
    (output / "meshioplusplus_fuzz_read_vtk").touch()
    with pytest.raises(ValueError, match="LFS pointer"):
        seeds.package(source, output)


@pytest.mark.parametrize("libraries", [False, True])
def test_campaign_library_opt_in(tmp_path, libraries):
    build, corpus, output = [tmp_path / name for name in ("build", "seeds", "out")]
    build.mkdir()
    corpus.mkdir()
    for name in ("meshioplusplus_fuzz_read", "meshioplusplus_fuzz_hdf5"):
        binary = build / name
        binary.write_text("#!/bin/sh\necho 'stat::number_of_executed_units: 1'\n")
        binary.chmod(0o755)
    env = dict(os.environ, FORMATS="vtk vtkhdf exodus vtx szplt xdmf elmer")
    if libraries:
        env.update(MIO_FUZZ_LIBRARY_READERS="1", FUZZ_BINARY="meshioplusplus_fuzz_hdf5")
    else:
        env["MIO_FUZZ_LIBRARY_READERS"] = "0"
        env["FUZZ_BINARY"] = "meshioplusplus_fuzz_read"
    subprocess.run(
        [
            "bash",
            str(REPO / "tools/fuzz/run_campaign.sh"),
            str(build),
            str(corpus),
            str(output),
            "1",
        ],
        env=env,
        check=True,
    )
    assert (output / "vtk/fuzz.log").exists()
    for fmt in ("vtkhdf", "exodus"):
        assert (output / fmt / "fuzz.log").exists() == libraries
    for fmt in ("vtx", "szplt", "xdmf", "elmer"):
        assert not (output / fmt).exists()


def test_native_connectivity_shape_probe(tmp_path):
    executable = os.environ.get("MIO_FUZZ_REPLAY")
    if not executable:
        pytest.skip("set MIO_FUZZ_REPLAY to test the native sanitizer harness")
    h5py = pytest.importorskip("h5py")
    import numpy as np

    path = tmp_path / "short-connectivity"
    with h5py.File(path, "w") as file:
        group = file.create_group("VTKHDF")
        group.attrs["Type"] = np.bytes_("UnstructuredGrid")
        group.attrs["Version"] = np.array([2, 0], dtype="int64")
        group["Points"] = np.array(
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        )
        for name, values in {
            "NumberOfPoints": [4],
            "NumberOfCells": [1],
            "NumberOfConnectivityIds": [4],
            "Connectivity": [0, 1, 2],
            "Offsets": [0, 4],
            "Types": [10],
        }.items():
            group[name] = np.array(values, dtype="int64")
    subprocess.run([executable, "-format=vtkhdf", str(path)], check=True, timeout=20)
