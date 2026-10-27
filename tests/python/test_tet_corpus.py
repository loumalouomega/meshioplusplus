"""tools/tet_corpus.py (roadmap 6.7.8): the sample drawn from the metadata,
the fetch's extraction and hash checks, one run end to end (with Netgen) and
the summary. No network: the dataset is a local archive and file:// URLs."""

import importlib.util
import io
import json
import pathlib
import tarfile

import pytest

import meshioplusplus as mio

try:
    from meshioplusplus import _core
except ImportError:  # pragma: no cover - the numpy-only install
    _core = None

ROOT = pathlib.Path(__file__).resolve().parents[2]
_spec = importlib.util.spec_from_file_location(
    "tet_corpus", ROOT / "tools" / "tet_corpus.py"
)
corpus = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(corpus)

needs_netgen = pytest.mark.skipif(
    not getattr(_core, "__has_netgen__", False), reason="needs the Netgen backend"
)

GEOMETRY = (
    "file_id,num_boundary_edges,edge_manifold,vertex_manifold,oriented\n"
    "1,0,1,1,1\n"
    "2,0,1,1,1\n"
    "3,12,1,1,1\n"
    "4,0,0,1,1\n"
    "5,0,1,1,1\n"
)
SUMMARY = (
    "ID,Link\n"
    "1,https://example.org/a.stl\n"
    "2,https://example.org/b.STL\n"
    "3,https://example.org/c.stl\n"
    "4,https://example.org/d.stl\n"
    "5,https://example.org/e.obj\n"
)


def test_the_sample_is_the_closed_manifold_stl_models_drawn_by_the_seed():
    import csv

    geometry = list(csv.DictReader(io.StringIO(GEOMETRY)))
    summary = list(csv.DictReader(io.StringIO(SUMMARY)))
    assert corpus.eligible(geometry, summary) == [1, 2]
    ids = list(range(100))
    a = corpus.draw(ids, 6708, 10)
    assert a == corpus.draw(list(reversed(ids)), 6708, 10) == sorted(a)
    assert len(set(a)) == 10 and a != corpus.draw(ids, 1, 10)
    assert corpus.draw([3, 1], 0, 10) == [1, 3]


def _dataset(tmp_path):
    """A stand-in for the dataset: its metadata and an archive of two STL
    models (a closed sphere and an open half of it) plus one left out."""
    base = tmp_path / "hf"
    (base / "metadata").mkdir(parents=True)
    (base / "metadata" / "geometry_data.csv").write_text(GEOMETRY)
    (base / "metadata" / "input_summary.csv").write_text(SUMMARY)
    sphere = mio.icosa_sphere(subdivisions=2)
    src = tmp_path / "raw"
    src.mkdir()
    mio.write(src / "1.stl", sphere)
    tri = sphere.cells[0].data
    keep = sphere.points[tri].mean(axis=1)[:, 2] > 0.0
    mio.write(src / "2.stl", mio.Mesh(sphere.points, [("triangle", tri[keep])]))
    mio.write(src / "3.stl", sphere)
    with tarfile.open(base / "Thingi10K.tar.gz", "w:gz") as tar:
        for name in ("1.stl", "2.stl", "3.stl"):
            tar.add(src / name, arcname=f"Thingi10K/raw_meshes/{name}")
    manifest = tmp_path / "manifest.json"
    manifest.write_text(
        json.dumps(
            {
                "Version": 1,
                "models": {},
                "sample": {"count": 5, "seed": 6708},
                "source": {
                    "archive": "Thingi10K.tar.gz",
                    "base": base.as_uri(),
                    "geometry": "metadata/geometry_data.csv",
                    "members": "Thingi10K/raw_meshes",
                    "summary": "metadata/input_summary.csv",
                },
            }
        )
    )
    return manifest, base


def test_fetch_extracts_the_sample_and_records_then_checks_its_hashes(tmp_path):
    manifest, base = _dataset(tmp_path)
    dest, record = tmp_path / "d", tmp_path / "hashes.json"
    argv = ["fetch", "--manifest", str(manifest), "--dest", str(dest)]
    assert corpus.main(argv + ["--record", str(record)]) == 0
    assert sorted(p.name for p in dest.glob("*.stl")) == ["1.stl", "2.stl"]
    assert json.loads((dest / "sample.json").read_text())["ids"] == [1, 2]
    hashes = json.loads(record.read_text())
    assert hashes == {
        "1": corpus.sha256(dest / "1.stl"),
        "2": corpus.sha256(dest / "2.stl"),
    }
    # A recorded hash that no longer matches removes the model and fails.
    doc = json.loads(manifest.read_text())
    doc["models"] = {"1": "0" * 64, "2": hashes["2"]}
    manifest.write_text(json.dumps(doc))
    assert corpus.main(argv) == 1
    assert sorted(p.name for p in dest.glob("*.stl")) == ["2.stl"]
    # A local copy of the archive is read in place of the URL.
    other = tmp_path / "e"
    doc["models"] = hashes
    manifest.write_text(json.dumps(doc))
    archive = str(base / "Thingi10K.tar.gz")
    assert (
        corpus.main(
            [
                "fetch",
                "--manifest",
                str(manifest),
                "--dest",
                str(other),
                "--archive",
                archive,
            ]
        )
        == 0
    )
    assert (other / "1.stl").exists()


@needs_netgen
def test_one_sample_end_to_end(tmp_path, capsys):
    manifest, _ = _dataset(tmp_path)
    dest, out = tmp_path / "d", tmp_path / "run.json"
    assert corpus.main(["fetch", "--manifest", str(manifest), "--dest", str(dest)]) == 0
    assert corpus.main(["run", "--dest", str(dest), "--out", str(out)]) == 0
    run = json.loads(out.read_text())
    rows = {r["id"]: r for r in run["models"]}
    assert rows[1]["status"] == "ok", rows[1]
    assert rows[1]["boundary_preserved"] is True and rows[1]["num_tetrahedra"] > 0
    assert 0.0 < rows[1]["min_dihedral_deg"] < rows[1]["max_dihedral_deg"] < 180.0
    assert rows[2]["status"] == "refused" and "open" in rows[2]["error"]
    assert corpus.main(["summary", str(out)]) == 0
    text = capsys.readouterr().out
    assert "| meshed | 1 of 1 (100.0%) |" in text
    assert "| accepted (closed, manifold, oriented as read) | 1 of 2 (50.0%) |" in text


def test_a_crash_fails_the_run_naming_the_model(tmp_path, monkeypatch, capsys):
    dest = tmp_path / "d"
    dest.mkdir()
    (dest / "sample.json").write_text(json.dumps({"ids": [7]}))
    (dest / "7.stl").write_text("solid x\nendsolid x\n")

    class Done:
        returncode = -11
        stdout = ""
        stderr = "Segmentation fault"

    monkeypatch.setattr(corpus.subprocess, "run", lambda *a, **k: Done())
    out = tmp_path / "run.json"
    assert corpus.main(["run", "--dest", str(dest), "--out", str(out)]) == 1
    assert "model 7 crashed or hung" in capsys.readouterr().err
    assert json.loads(out.read_text())["models"][0]["status"] == "crashed"
