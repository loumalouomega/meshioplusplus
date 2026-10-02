"""Caller-buffer reports run once, stay owned, and match the Python surface."""

import ctypes
import json
import os
from pathlib import Path

import pytest

import meshioplusplus as mp

from . import helpers


def _library():
    path = Path(
        os.environ.get(
            "MESHIOPLUSPLUS_LIB",
            Path(__file__).resolve().parents[2]
            / "build/cpp-release/libmeshioplusplus.so",
        )
    )
    if not path.is_file():
        pytest.skip("requires the native C library")
    lib = ctypes.CDLL(str(path))
    if not hasattr(lib, "mio_pipeline_report_json"):
        pytest.skip("requires the rebuilt report API")
    for name in (
        "mio_pipeline_run_json_report",
        "mio_pipeline_run_file_report",
        "mio_sequence_pipeline_run_json_report",
        "mio_sequence_pipeline_run_file_report",
    ):
        fn = getattr(lib, name)
        fn.argtypes = [ctypes.c_char_p]
        fn.restype = ctypes.c_void_p
    lib.mio_pipeline_report_json.argtypes = [
        ctypes.c_void_p,
        ctypes.c_void_p,
        ctypes.c_int64,
    ]
    lib.mio_pipeline_report_json.restype = ctypes.c_int64
    lib.mio_pipeline_report_free.argtypes = [ctypes.c_void_p]
    lib.mio_last_error.restype = ctypes.c_char_p
    return lib


@pytest.mark.parametrize("file_input", [False, True])
def test_c_pipeline_report_matches_python_and_survives_input_removal(
    tmp_path, file_input
):
    lib = _library()
    if not lib.mio_pipeline_has_json():
        pytest.skip("requires JSON support")
    source = tmp_path / "in.vtu"
    target = tmp_path / "out.vtu"
    mp.write(source, helpers.tet_mesh, compression=None)
    settings = {
        "Input": {"Path": str(source)},
        "Output": {"Path": str(target), "Codec": "none"},
        "Operations": [{"Op": "Quality"}, {"Op": "Transform", "Translate": [1, 0, 0]}],
    }
    expected = mp.run_pipeline(settings)
    text = json.dumps(settings)
    if file_input:
        path = tmp_path / "settings.json"
        path.write_text(text)
        handle = lib.mio_pipeline_run_file_report(os.fsencode(path))
    else:
        handle = lib.mio_pipeline_run_json_report(text.encode())
    assert handle, lib.mio_last_error()
    try:
        source.unlink()  # Reading the report must not run the pipeline again.
        n = lib.mio_pipeline_report_json(handle, None, 0)
        assert n > 0
        small = ctypes.create_string_buffer(2)
        assert lib.mio_pipeline_report_json(handle, small, 2) == n
        assert small.value == b"{"
        out = ctypes.create_string_buffer(n + 1)
        assert lib.mio_pipeline_report_json(handle, out, len(out)) == n
        assert json.loads(out.value) == expected
        assert lib.mio_pipeline_report_json(handle, None, -1) == -1
    finally:
        lib.mio_pipeline_report_free(handle)
    lib.mio_pipeline_report_free(None)
    assert lib.mio_pipeline_report_json(None, None, 0) == -1


def test_c_report_error_guards_and_disabled_capability():
    lib = _library()
    for name in (
        "mio_pipeline_run_json_report",
        "mio_pipeline_run_file_report",
        "mio_sequence_pipeline_run_json_report",
        "mio_sequence_pipeline_run_file_report",
    ):
        assert getattr(lib, name)(None) is None
        assert b"NULL" in lib.mio_last_error()
    assert lib.mio_pipeline_run_json_report(b"{}") is None
    if not lib.mio_pipeline_has_json():
        assert b"MESHIOPLUSPLUS_WITH_JSON" in lib.mio_last_error()
    else:
        assert b"Input" in lib.mio_last_error()


@pytest.mark.parametrize("file_input", [False, True])
def test_c_sequence_report_shape(tmp_path, file_input):
    lib = _library()
    if not lib.mio_pipeline_has_json():
        pytest.skip("requires JSON support")
    source = tmp_path / "in.vtu"
    mp.write(source, helpers.tet_mesh, compression=None)
    target = str(tmp_path / "out_{step}.vtu")
    settings = {
        "Input": {"Paths": [str(source)]},
        "Output": {"Path": target, "Codec": "none"},
        "Operations": [{"Op": "Quality"}],
    }
    text = json.dumps(settings)
    if file_input:
        path = tmp_path / "settings.json"
        path.write_text(text)
        handle = lib.mio_sequence_pipeline_run_file_report(os.fsencode(path))
    else:
        handle = lib.mio_sequence_pipeline_run_json_report(text.encode())
    assert handle, lib.mio_last_error()
    try:
        n = lib.mio_pipeline_report_json(handle, None, 0)
        buf = ctypes.create_string_buffer(n + 1)
        assert lib.mio_pipeline_report_json(handle, buf, len(buf)) == n
        report = json.loads(buf.value)
        assert report["steps"][0]["op"] == "Quality"
        assert len(report["warnings"]) == 1
        assert "using the integer index" in report["warnings"][0]
        assert list(tmp_path.glob("out_*.vtu"))
    finally:
        lib.mio_pipeline_report_free(handle)
