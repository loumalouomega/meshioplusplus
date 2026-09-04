"""The HTTP front-end of the MCP server (``meshioplusplus-mcp --http``).

Gated on the ``[dashboard]`` extra (Starlette/uvicorn, plus httpx for
Starlette's test client): the pure tool layer these routes dispatch to is
covered by ``test_mcp.py`` in the default matrix. Every test here drives the
real ASGI app through ``TestClient`` -- no socket, no uvicorn.
"""

from __future__ import annotations

import pytest

pytest.importorskip("mcp")
pytest.importorskip("starlette")
pytest.importorskip("httpx")

import numpy as np  # noqa: E402
from starlette.testclient import TestClient  # noqa: E402

import meshioplusplus  # noqa: E402
from meshioplusplus.mcp import _tools, create_server, has_dashboard  # noqa: E402
from meshioplusplus.mcp._http import (  # noqa: E402
    DEFAULT_ALLOWED_ORIGINS,
    build_app,
    default_runs_dir,
    new_token,
    origin_regex,
)


def _mixed_mesh():
    points = np.array(
        [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1]], dtype=float
    )
    cells = [
        ("triangle", np.array([[0, 1, 2], [0, 2, 3]])),
        ("tetra", np.array([[0, 1, 2, 4]])),
    ]
    return meshioplusplus.Mesh(
        points, cells, point_data={"t": np.linspace(0.0, 1.0, 5)}
    )


TOKEN = "s3cret-token"
PAGES = "https://loumalouomega.github.io"


@pytest.fixture(autouse=True)
def _reset_root():
    yield
    _tools.set_root(None)


def _client(root, token=TOKEN, **kwargs):
    app = build_app(
        create_server(root=str(root)), token=token, root=str(root), **kwargs
    )
    return TestClient(app)


def _auth(token=TOKEN):
    return {"Authorization": f"Bearer {token}"}


def test_has_dashboard_here():
    assert has_dashboard()
    assert len(new_token()) >= 24


def test_health_requires_the_token_and_describes_the_server(tmp_path):
    with _client(tmp_path) as client:
        denied = client.get("/api/health")
        assert denied.status_code == 401
        assert denied.json()["error_type"] == "PermissionError"
        wrong = client.get("/api/health", headers=_auth("nope"))
        assert wrong.status_code == 401
        ok = client.get("/api/health", headers=_auth())
        assert ok.status_code == 200
        body = ok.json()
        assert body["version"] == meshioplusplus.__version__
        assert body["root"] == str(tmp_path)
        assert body["runs_dir"] == default_runs_dir(str(tmp_path))
        assert body["auth"] == "token"
        assert body["mcp"] == "/mcp"
        assert body["transport"] in ("streamable-http", "sse")
        assert "dataset_health" in body["tools"] and "info" in body["tools"]
        assert body["tools"] == list(_tools.TOOL_REGISTRY)


def test_no_token_mode(tmp_path):
    with _client(tmp_path, token=None) as client:
        body = client.get("/api/health").json()
        assert body["auth"] == "none"


def test_tools_dispatch_and_error_shapes(tmp_path):
    mesh_file = tmp_path / "in.vtu"
    meshioplusplus.write(str(mesh_file), _mixed_mesh())
    with _client(tmp_path) as client:
        ok = client.post(
            "/api/tools/info", json={"input_path": "in.vtu"}, headers=_auth()
        )
        assert ok.status_code == 200
        assert ok.json()["num_points"] == 5
        # a tool failure is a 200 payload, the MCP rule
        failed = client.post(
            "/api/tools/info", json={"input_path": "nope.vtu"}, headers=_auth()
        )
        assert failed.status_code == 200
        assert failed.json()["error_type"] == "ValueError"
        # paths inside the body go through the tools' own sandbox
        outside = client.post(
            "/api/tools/info", json={"input_path": "../outside.vtu"}, headers=_auth()
        )
        assert "outside the configured root" in outside.json()["error"]
        unknown = client.post("/api/tools/nope", json={}, headers=_auth())
        assert unknown.status_code == 404
        assert unknown.json()["error_type"] == "KeyError"
        bad = client.post(
            "/api/tools/info",
            content=b"[",
            headers={**_auth(), "Content-Type": "application/json"},
        )
        assert bad.status_code == 400
        not_object = client.post("/api/tools/info", json=[1, 2], headers=_auth())
        assert not_object.status_code == 400
        assert not_object.json()["error_type"] == "TypeError"
        empty = client.post("/api/tools/formats", headers=_auth())
        assert empty.status_code == 200 and empty.json()["readable"]
        unauth = client.post("/api/tools/info", json={"input_path": "in.vtu"})
        assert unauth.status_code == 401


def test_files_are_sandboxed_downloads(tmp_path):
    payload = b"checkpoint bytes"
    (tmp_path / "runs").mkdir()
    (tmp_path / "runs" / "best.mdlus").write_bytes(payload)
    (tmp_path.parent / "secret.txt").write_bytes(b"no")
    with _client(tmp_path) as client:
        ok = client.get(
            "/api/files", params={"path": "runs/best.mdlus"}, headers=_auth()
        )
        assert ok.status_code == 200 and ok.content == payload
        assert "best.mdlus" in ok.headers["content-disposition"]
        # a download link cannot set a header: ?token= is accepted here only
        via_query = client.get(
            "/api/files", params={"path": "runs/best.mdlus", "token": TOKEN}
        )
        assert via_query.status_code == 200
        via_query_elsewhere = client.get("/api/health", params={"token": TOKEN})
        assert via_query_elsewhere.status_code == 401
        outside = client.get(
            "/api/files", params={"path": "../secret.txt"}, headers=_auth()
        )
        assert outside.status_code == 404
        assert "outside the configured root" in outside.json()["error"]
        missing = client.get("/api/files", headers=_auth())
        assert missing.status_code == 400
        directory = client.get("/api/files", params={"path": "runs"}, headers=_auth())
        assert directory.status_code == 404


def _preflight(client, origin, extra=None):
    return client.options(
        "/api/tools/info",
        headers={
            "Origin": origin,
            "Access-Control-Request-Method": "POST",
            "Access-Control-Request-Headers": "authorization,content-type",
            **(extra or {}),
        },
    )


def test_cors_admits_loopback_and_the_docs_site_only(tmp_path):
    with _client(tmp_path) as client:
        for origin in ("http://localhost:5173", "http://127.0.0.1:4173", PAGES):
            response = _preflight(client, origin)
            assert response.status_code == 200, origin
            assert response.headers["access-control-allow-origin"] == origin
            assert (
                "authorization"
                in response.headers["access-control-allow-headers"].lower()
            )
        evil = _preflight(client, "https://evil.example")
        assert "access-control-allow-origin" not in evil.headers
        # an extra origin can be admitted explicitly
    with _client(
        tmp_path, allowed_origins=[*DEFAULT_ALLOWED_ORIGINS, "https://intranet.example"]
    ) as client:
        response = _preflight(client, "https://intranet.example")
        assert (
            response.headers["access-control-allow-origin"]
            == "https://intranet.example"
        )
        # a 401 still carries CORS headers, so the browser sees a 401, not an opaque error
        denied = client.get("/api/health", headers={"Origin": PAGES})
        assert denied.status_code == 401
        assert denied.headers["access-control-allow-origin"] == PAGES


def test_private_network_access_preflight_is_answered(tmp_path):
    with _client(tmp_path) as client:
        response = _preflight(
            client, PAGES, {"Access-Control-Request-Private-Network": "true"}
        )
        # Starlette >= 0.46 rejects a PNA preflight outright unless allowed --
        # the status is the assertion that matters, not just the header.
        assert response.status_code == 200, response.text
        assert response.headers["access-control-allow-private-network"] == "true"
        plain = _preflight(client, PAGES)
        assert "access-control-allow-private-network" not in plain.headers


def test_origin_regex():
    import re

    pattern = re.compile(origin_regex(["https://a.example/"]))
    assert pattern.match("http://localhost:5173")
    assert pattern.match("http://127.0.0.1")
    assert pattern.match("https://a.example")
    assert not pattern.match("https://a.example.evil")
    assert not pattern.match("http://localhost.evil:80")


def test_main_prints_the_named_error_without_the_extra(monkeypatch, capsys):
    import meshioplusplus.mcp as mmcp
    from meshioplusplus.mcp import _server

    monkeypatch.setattr(mmcp, "has_dashboard", lambda: False)
    assert _server.main(["--http"]) == 1
    assert "meshioplusplus[dashboard]" in capsys.readouterr().err
