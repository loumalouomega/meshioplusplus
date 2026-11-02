"""Meshes read from URLs through the optional fsspec (roadmap 7.1.2):
``read``, ``read_metadata``, sequences and dataset manifests, over fsspec's
in-memory filesystem, ``file://`` and a local HTTP server -- no network."""

import contextlib
import functools
import gzip
import http.server
import os
import threading

import numpy as np
import pytest

import meshioplusplus as mp
from meshioplusplus import _remote

from . import helpers

fsspec = pytest.importorskip("fsspec")


def _vtu_bytes(tmp_path, mesh=None, name="cube.vtu"):
    path = tmp_path / name
    mp.write(str(path), mesh if mesh is not None else helpers.tet_mesh)
    return path.read_bytes()


@pytest.fixture
def memfs():
    fs = fsspec.filesystem("memory")
    yield fs
    for path in fs.find("/"):
        fs.rm(path)


def _same(a, b):
    np.testing.assert_allclose(a.points, b.points, atol=1e-12)
    assert [c.type for c in a.cells] == [c.type for c in b.cells]
    for x, y in zip(a.cells, b.cells):
        np.testing.assert_array_equal(x.data, y.data)


def test_is_url():
    assert _remote.is_url("https://example.org/cube.vtu")
    assert _remote.is_url("s3://bucket/run/step.vtu")
    assert _remote.is_url("simplecache::https://example.org/cube.vtu")
    assert _remote.is_url("file:///tmp/cube.vtu")
    assert not _remote.is_url("cube.vtu")
    assert not _remote.is_url("C:\\data\\cube.vtu")
    assert not _remote.is_url("C://data/cube.vtu")
    assert not _remote.is_url(os.path.join("a", "b.vtu"))
    assert _remote.url_name("https://h/x/cube.vtu?sig=abc#frag") == "cube.vtu"
    assert _remote.url_name("https://h/a%20b.vtu") == "a b.vtu"
    assert _remote.url_name("https://h/") == "download"


def test_read_from_memory_and_file_urls(tmp_path, memfs):
    data = _vtu_bytes(tmp_path)
    memfs.pipe("/run/cube.vtu", data)
    local = mp.read(str(tmp_path / "cube.vtu"))
    _same(mp.read("memory://run/cube.vtu"), local)
    _same(mp.read((tmp_path / "cube.vtu").as_uri()), local)


def test_compressed_url(tmp_path, memfs):
    memfs.pipe("/run/cube.vtu.gz", gzip.compress(_vtu_bytes(tmp_path)))
    _same(mp.read("memory://run/cube.vtu.gz"), mp.read(str(tmp_path / "cube.vtu")))


def test_read_metadata_from_url(tmp_path, memfs):
    memfs.pipe("/run/cube.vtu", _vtu_bytes(tmp_path))
    meta = mp.read_metadata("memory://run/cube.vtu")
    assert meta["num_points"] == len(helpers.tet_mesh.points)
    assert meta["format"] == "vtu"


def test_missing_url_and_multi_file_format(memfs):
    with pytest.raises(mp.ReadError, match="not found"):
        mp.read("memory://nothing/here.vtu")
    memfs.pipe("/run/case.pvtu", b"<VTKFile/>")
    with pytest.raises(mp.ReadError, match="pvtu.*several files"):
        mp.read("memory://run/case.pvtu")
    with pytest.raises(mp.ReadError, match="tetgen.*several files"):
        mp.read("memory://run/case.pvtu", file_format="tetgen")


def test_unknown_protocol_names_the_package():
    with pytest.raises(mp.ReadError, match="cannot read|not found|protocol"):
        mp.read("nosuchprotocol://bucket/cube.vtu")


def test_sequence_over_a_url_pattern(tmp_path, memfs):
    for k in (2, 10, 1):
        mesh = mp.Mesh(helpers.tet_mesh.points + k, helpers.tet_mesh.cells)
        memfs.pipe(f"/series/step_{k}.vtu", _vtu_bytes(tmp_path, mesh, f"s{k}.vtu"))
    entries = mp.sequence_entries("memory://series/step_*.vtu")
    assert [os.path.basename(e["path"]) for e in entries] == [
        "step_1.vtu",
        "step_2.vtu",
        "step_10.vtu",
    ]
    assert [e["time"] for e in entries] == [1.0, 2.0, 10.0]
    steps = list(mp.read_sequence("memory://series/step_*.vtu"))
    assert [t for t, _ in steps] == [1.0, 2.0, 10.0]
    np.testing.assert_allclose(steps[2][1].points, helpers.tet_mesh.points + 10)
    with pytest.raises(mp.ReadError, match="matched no files"):
        mp.sequence_entries("memory://series/none_*.vtu")
    with pytest.raises(mp.ReadError, match="not found"):
        mp.sequence_entries(["memory://series/step_1.vtu", "memory://series/gone.vtu"])


def test_dataset_manifest_keeps_url_sources(tmp_path):
    from meshioplusplus._dataset import _resolve_source_path

    url = "s3://bucket/cases/run_*.vtu"
    assert _resolve_source_path(url, str(tmp_path)) == url
    assert _resolve_source_path("cases/a.vtu", "/base").endswith(
        os.path.join("base", "cases", "a.vtu")
    )


@contextlib.contextmanager
def _http_server(directory):
    handler = functools.partial(
        http.server.SimpleHTTPRequestHandler, directory=str(directory)
    )
    handler.log_message = lambda *a, **k: None
    server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        yield f"http://127.0.0.1:{server.server_address[1]}"
    finally:
        server.shutdown()
        server.server_close()


def test_read_over_http(tmp_path):
    pytest.importorskip("aiohttp")
    _vtu_bytes(tmp_path)
    with _http_server(tmp_path) as base:
        _same(mp.read(f"{base}/cube.vtu"), mp.read(str(tmp_path / "cube.vtu")))
        with pytest.raises(mp.ReadError, match="not found"):
            mp.read(f"{base}/missing.vtu")
