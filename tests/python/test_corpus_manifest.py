"""The fixture manifest lists every reference fixture and nothing else."""

import hashlib
import importlib.util
import json
import pathlib

REPO = pathlib.Path(__file__).resolve().parents[2]
_spec = importlib.util.spec_from_file_location(
    "gen_corpus_manifest", REPO / "tools" / "gen_corpus_manifest.py"
)
gen = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(gen)


def _rows():
    return json.loads(gen.MANIFEST.read_text(encoding="utf-8"))["fixtures"]


def test_every_fixture_is_listed():
    listed = {row["path"] for row in _rows()}
    on_disk = set(gen.fixture_paths())
    assert (
        on_disk - listed == set()
    ), "fixture without a manifest row: run tools/gen_corpus_manifest.py"


def test_every_row_exists():
    missing = [row["path"] for row in _rows() if not (REPO / row["path"]).is_file()]
    assert (
        missing == []
    ), "manifest row without a file: run tools/gen_corpus_manifest.py"


def test_the_manifest_is_fresh():
    assert gen.MANIFEST.read_text(encoding="utf-8") == gen.render(gen.build())


def test_fetched_fixtures_match_their_hash():
    for row in _rows():
        data = (REPO / row["path"]).read_bytes()
        if data.startswith(b"version https://git-lfs.github.com/spec/v1"):
            continue  # not fetched; the manifest holds the pointer's oid
        assert hashlib.sha256(data).hexdigest() == row["sha256"], row["path"]


def test_documented_directories_point_at_an_existing_readme():
    for row in _rows():
        if row["readme"] is not None:
            assert (REPO / row["readme"]).is_file()
