"""Tests for corpus packaging and the opt-in library campaign."""

import importlib.util
import os
import pathlib
import shutil
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


def test_packages_generated_and_regression_seeds(tmp_path):
    source, generated, output = [
        tmp_path / name for name in ("regressions", "generated", "out")
    ]
    for root in (source, generated):
        (root / "vtu").mkdir(parents=True)
    output.mkdir()
    (source / "vtu/regression").write_bytes(b"malformed input")
    (generated / "vtu/positive.vtu").write_bytes(b"valid input")
    (generated / "vtu/empty").touch()
    (generated / "vtu/subdirectory").mkdir()
    (output / "meshioplusplus_fuzz_read_vtu").touch()
    (output / "meshioplusplus_fuzz_read_vtu.dict").touch()
    (output / "meshioplusplus_fuzz_read_off").touch()
    seeds.package(source, output, generated=generated)
    with zipfile.ZipFile(
        output / "meshioplusplus_fuzz_read_vtu_seed_corpus.zip"
    ) as archive:
        assert archive.namelist() == ["regression", "positive.vtu"]
        assert archive.read("positive.vtu") == b"valid input"
    assert not (output / "meshioplusplus_fuzz_read_off_seed_corpus.zip").exists()


@pytest.mark.parametrize("list_fails", [False, True])
def test_oss_fuzz_build_publishes_registry_targets(tmp_path, list_fails):
    tools, work, output = [tmp_path / name for name in ("bin", "work", "out")]
    for directory in (tools, work, output):
        directory.mkdir()
    cmake = tools / "cmake"
    cmake.write_text("#!/bin/sh\nexit 0\n")
    cmake.chmod(0o755)
    build = work / "meshioplusplus-build"
    build.mkdir()
    scripts = {
        "meshioplusplus_fuzz_read": "#!/bin/sh\nexit 0\n",
        "meshioplusplus_fuzz_replay": (
            "#!/bin/sh\nexit 1\n"
            if list_fails
            else "#!/bin/sh\nprintf 'vtu\\nno_dictionary\\nelmer\\ncgns\\n'\n"
        ),
        "meshioplusplus_fuzz_seeds": (
            '#!/bin/sh\nmkdir -p "$1/vtu"\nprintf positive > "$1/vtu/seed.vtu"\n'
        ),
    }
    for name, text in scripts.items():
        executable = build / name
        executable.write_text(text)
        executable.chmod(0o755)
    source = tmp_path / "src"
    fuzz = source / "meshioplusplus/tools/fuzz"
    (fuzz / "dicts").mkdir(parents=True)
    (fuzz / "dicts/vtu.dict").write_text('"VTKFile"\n')
    (fuzz / "not_fuzzed.txt").write_text("elmer\ncgns\n")
    (fuzz / "oss-fuzz").mkdir()
    shutil.copyfile(
        REPO / "tools/fuzz/oss-fuzz/package_seeds.py",
        fuzz / "oss-fuzz/package_seeds.py",
    )
    env = dict(
        os.environ,
        PATH=f"{tools}{os.pathsep}{os.environ['PATH']}",
        SRC=str(source),
        WORK=str(work),
        OUT=str(output),
        LIB_FUZZING_ENGINE="-fsanitize=fuzzer",
    )
    result = subprocess.run(
        ["bash", "-eu", str(REPO / "tools/fuzz/oss-fuzz/build.sh")],
        env=env,
        capture_output=True,
        text=True,
    )
    if list_fails:
        assert result.returncode != 0
        assert not list(output.iterdir())
    else:
        assert result.returncode == 0, result.stderr
        for fmt in ("vtu", "no_dictionary"):
            assert (output / f"meshioplusplus_fuzz_read_{fmt}").exists()
        for fmt in ("elmer", "cgns"):
            assert not (output / f"meshioplusplus_fuzz_read_{fmt}").exists()
        assert (output / "meshioplusplus_fuzz_read_vtu.dict").exists()
        assert not (output / "meshioplusplus_fuzz_read_no_dictionary.dict").exists()
        with zipfile.ZipFile(
            output / "meshioplusplus_fuzz_read_vtu_seed_corpus.zip"
        ) as archive:
            assert archive.read("seed.vtu") == b"positive"


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
    for fmt in ("vtx", "szplt", "xdmf"):
        assert (output / fmt / "fuzz.log").exists()
    assert not (output / "elmer").exists()


def test_bundle_pack_round_trip():
    import struct

    def pack(entries):
        out = bytearray(b"MIOB\x01") + struct.pack("<H", len(entries))
        for name, data in entries:
            name_b = name.encode()
            out += struct.pack("<H", len(name_b)) + name_b
            out += struct.pack("<I", len(data)) + data
        return bytes(out)

    blob = pack([("input.xdmf", b"<Xdmf/>"), ("data.h5", b"\x89HDF\r\n\x1a\n")])
    assert blob.startswith(b"MIOB\x01")
    assert blob[5:7] == struct.pack("<H", 2)


def test_seed_corpus_packs_xdmf_bundles(tmp_path):
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "seed_corpus", REPO / "tools/fuzz/seed_corpus.py"
    )
    seed_corpus = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(seed_corpus)
    out = tmp_path / "seeds"
    out.mkdir()
    assert seed_corpus._xdmf_bundle_seeds(out, 262144) >= 2
    assert list(out.glob("bundle-*.miob"))
    assert list(out.glob("bundle-regions-*.miob"))


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
