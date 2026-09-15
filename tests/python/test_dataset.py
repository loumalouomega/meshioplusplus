"""Tests for dataset manifests (`_dataset.py`).

Everything here is pure stdlib + numpy and runs in the default CI matrix --
the manifest deliberately has no optional dependency. File-backed cases write
tiny real meshes into tmp_path via the ordinary `write`.
"""

from __future__ import annotations

import json
import os

import numpy as np
import pytest

import meshioplusplus
from meshioplusplus import DatasetEntry, DatasetManifest
from meshioplusplus._dataset import check_pairing
from meshioplusplus._mesh import Mesh


def _tiny_mesh(offset=0.0):
    points = (
        np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]])
        + offset
    )
    return Mesh(
        points,
        [("triangle", np.array([[0, 1, 2], [0, 2, 3]], dtype=np.int64))],
        point_data={"T": np.arange(4, dtype=np.float64) + offset},
    )


def _write_cases(directory, n=3):
    directory.mkdir(exist_ok=True)
    paths = []
    for i in range(n):
        path = directory / f"case_{i}.vtu"
        meshioplusplus.write(str(path), _tiny_mesh(float(i)))
        paths.append(str(path))
    return paths


# --------------------------------------------------------------------------- #
# load: the tri-modal rule and strict validation                              #
# --------------------------------------------------------------------------- #
def _doc(**overrides):
    doc = {
        "Version": 1,
        "Name": "campaign",
        "Entries": [
            {
                "Id": "a",
                "Source": {"Pattern": "cases/out_*.vtu"},
                "Split": "train",
                "Tags": ["re100"],
                "Group": "cyl/laminar",
                "Notes": "n",
                "Metadata": {"Re": 100},
            },
            {"Id": "b", "Source": {"Path": "cases/single.vtu"}},
        ],
    }
    doc.update(overrides)
    return doc


def test_load_accepts_dict_json_text_and_path(tmp_path):
    doc = _doc()
    from_dict = DatasetManifest.load(doc)
    from_text = DatasetManifest.load(json.dumps(doc))
    path = tmp_path / "m.json"
    path.write_text(json.dumps(doc), encoding="utf-8")
    from_path = DatasetManifest.load(path)
    for m in (from_dict, from_text, from_path):
        assert m.name == "campaign"
        assert m.ids() == ["a", "b"]
        assert m["a"].split == "train"
        assert m["a"].tags == ("re100",)
        assert m["a"].metadata == {"Re": 100}
    # only the file-backed load records a base_dir
    assert from_dict.base_dir is None
    assert from_path.base_dir == str(tmp_path)
    assert from_path["a"].base_dir == str(tmp_path)


@pytest.mark.parametrize(
    "mutate, needle",
    [
        (lambda d: d.update(Bogus=1), "unknown key 'Bogus'"),
        (lambda d: d.update(Version=2), "unsupported Version 2"),
        (lambda d: d["Entries"][0].pop("Id"), "Entries[0].Id"),
        (lambda d: d["Entries"][0].update(Extra=1), "unknown key 'Extra'"),
        (lambda d: d["Entries"][0].update(Source={}), "exactly one of"),
        (
            lambda d: d["Entries"][0].update(Source={"Pattern": "x", "Path": "y"}),
            "exactly one of",
        ),
        # A Target failure must name Target, not Source -- the validator is
        # shared, so the block name has to be threaded through it.
        (
            lambda d: d["Entries"][0].update(Target={"Pattern": "x", "Zz": 1}),
            "unknown key 'Zz' in Entries[0].Target",
        ),
        (
            lambda d: d["Entries"][0].update(Target={}),
            "Entries[0].Target needs exactly one of",
        ),
        (
            lambda d: d["Entries"][0].update(Source={"Pattern": "x", "Zz": 1}),
            "unknown key 'Zz'",
        ),
        (
            lambda d: d["Entries"][0].update(
                Source={"Pattern": "x", "TimeFrom": "bogus"}
            ),
            "TimeFrom",
        ),
        (
            lambda d: d["Entries"][0].update(Source={"Pattern": "x", "Sort": True}),
            "Sort applies only with Paths",
        ),
        (lambda d: d["Entries"][0].update(Tags=[1]), "Tags"),
        (lambda d: d["Entries"][1].update(Id="a"), "duplicate entry id 'a'"),
    ],
)
def test_strict_validation_names_the_offender(mutate, needle):
    doc = _doc()
    mutate(doc)
    with pytest.raises(ValueError, match="meshio\\+\\+: dataset:") as excinfo:
        DatasetManifest.load(doc)
    assert needle in str(excinfo.value)


def test_round_trip_is_exact_and_stable(tmp_path):
    doc = _doc()
    manifest = DatasetManifest.load(doc)
    assert manifest.to_dict() == doc
    path = tmp_path / "m.json"
    manifest.save(path)
    first = path.read_text(encoding="utf-8")
    assert first.endswith("\n")
    DatasetManifest.load(path).save(path)
    assert path.read_text(encoding="utf-8") == first  # byte-stable
    assert DatasetManifest.load(path).to_dict() == doc


def test_to_dict_omits_empty_optional_keys():
    manifest = DatasetManifest()
    manifest.add({"Path": "a.vtu"}, validate_source=False)
    entry_doc = manifest.to_dict()["Entries"][0]
    assert set(entry_doc) == {"Id", "Source"}


# --------------------------------------------------------------------------- #
# source resolution: relative paths, TimeSeries, plans                        #
# --------------------------------------------------------------------------- #
def test_relative_sources_resolve_against_the_manifest_dir(tmp_path):
    _write_cases(tmp_path / "cases")
    doc = {
        "Version": 1,
        "Entries": [{"Id": "c", "Source": {"Pattern": "cases/case_*.vtu"}}],
    }
    path = tmp_path / "sub" / ".." / "m.json"
    os.makedirs(os.path.dirname(path), exist_ok=True)
    real = tmp_path / "m.json"
    real.write_text(json.dumps(doc), encoding="utf-8")
    manifest = DatasetManifest.load(real)
    plan = manifest["c"].entries()
    assert len(plan) == 3
    assert all(os.path.isabs(e["path"]) for e in plan)
    assert plan[0]["path"].endswith(os.path.join("cases", "case_0.vtu"))

    series = manifest["c"].time_series()
    assert len(series) == 3
    t0, mesh0 = series[0]
    assert mesh0.points.shape == (4, 3)


def test_absolute_and_cwd_relative_sources(tmp_path, monkeypatch):
    paths = _write_cases(tmp_path / "cases")
    # absolute paths ignore base_dir entirely
    manifest = DatasetManifest(base_dir="/nonexistent")
    entry = manifest.add({"Paths": paths})
    assert len(entry.entries()) == 3
    # a dict-loaded manifest (no base_dir) resolves against the CWD
    monkeypatch.chdir(tmp_path)
    m2 = DatasetManifest.load(
        {"Version": 1, "Entries": [{"Id": "x", "Source": {"Path": "cases/case_1.vtu"}}]}
    )
    assert len(m2["x"].entries()) == 1


def test_source_plan_kwargs_reach_the_sequence_layer(tmp_path):
    paths = _write_cases(tmp_path / "cases")
    manifest = DatasetManifest()
    entry = manifest.add(
        {
            "Paths": [paths[2], paths[0], paths[1]],
            "Times": [0.5, 1.5, 2.5],
            "Sort": True,
        },
        id="c",
    )
    plan = entry.entries()
    # Sort=True natural-sorts the stated list; Times apply after ordering.
    assert [os.path.basename(e["path"]) for e in plan] == [
        "case_0.vtu",
        "case_1.vtu",
        "case_2.vtu",
    ]
    assert [e["time"] for e in plan] == [0.5, 1.5, 2.5]
    assert all(e["time_source"] == "explicit" for e in plan)
    series = entry.time_series()
    assert series.times == [0.5, 1.5, 2.5]


def test_add_validate_source_fails_early_on_an_empty_glob(tmp_path):
    manifest = DatasetManifest(base_dir=str(tmp_path))
    with pytest.raises(meshioplusplus.ReadError):
        manifest.add({"Pattern": "missing_*.vtu"})
    assert len(manifest) == 0  # nothing half-added
    entry = manifest.add({"Pattern": "missing_*.vtu"}, validate_source=False)
    assert entry.id == "missing"


# --------------------------------------------------------------------------- #
# add / remove / id derivation                                                #
# --------------------------------------------------------------------------- #
def test_add_shapes_and_id_derivation():
    manifest = DatasetManifest()
    a = manifest.add("out_*.vtu", validate_source=False)
    assert a.source == {"Pattern": "out_*.vtu"} and a.id == "out"
    b = manifest.add("single.vtu", validate_source=False)
    assert b.source == {"Path": "single.vtu"} and b.id == "single"
    c = manifest.add(["p0.vtu", "p1.vtu"], validate_source=False)
    assert c.source == {"Paths": ["p0.vtu", "p1.vtu"]} and c.id == "p0"
    with pytest.raises(ValueError, match="already exists"):
        manifest.add("out_9*.vtu", id="out", validate_source=False)
    with pytest.raises(ValueError, match="already exists"):
        manifest.add("out_*.vtu", validate_source=False)  # derived collision
    manifest.remove("out")
    assert manifest.ids() == ["single", "p0"]
    with pytest.raises(KeyError):
        manifest["out"]


# --------------------------------------------------------------------------- #
# curation: split / assign_splits / tag / annotate                            #
# --------------------------------------------------------------------------- #
def _many(n=10, groups=False):
    manifest = DatasetManifest()
    for i in range(n):
        manifest.add(
            f"case_{i}.vtu",
            id=f"c{i}",
            group=f"g{i // 2}" if groups else None,
            validate_source=False,
        )
    return manifest


def test_set_split_and_splits():
    manifest = _many(4)
    manifest.set_split(["c0", "c1"], "train")
    manifest.set_split("c2", "test")
    assert manifest.splits() == {"train": 2, "test": 1, None: 1}
    manifest.set_split("all", None)
    assert manifest.splits() == {None: 4}
    with pytest.raises(KeyError):
        manifest.set_split(["c0", "nope"], "train")


def test_assign_splits_is_deterministic_and_exhaustive():
    fractions = {"train": 0.8, "valid": 0.1, "test": 0.1}
    a, b = _many(20), _many(20)
    a.assign_splits(fractions, seed=7)
    b.assign_splits(fractions, seed=7)
    assert [e.split for e in a] == [e.split for e in b]
    counts = a.splits()
    assert None not in counts and sum(counts.values()) == 20
    assert counts["train"] == 16 and counts["valid"] == 2 and counts["test"] == 2
    c = _many(20)
    c.assign_splits(fractions, seed=8)
    assert [e.split for e in a] != [e.split for e in c]  # the seed matters


def test_assign_splits_by_group_keeps_groups_together():
    manifest = _many(10, groups=True)
    manifest.assign_splits({"train": 0.6, "test": 0.4}, seed=0, by_group=True)
    by_group = {}
    for entry in manifest:
        by_group.setdefault(entry.group, set()).add(entry.split)
    assert all(len(s) == 1 for s in by_group.values())
    assert sum(manifest.splits().values()) == 10


def test_assign_splits_validates_fractions():
    manifest = _many(4)
    with pytest.raises(ValueError, match="sum to 1"):
        manifest.assign_splits({"train": 0.5, "test": 0.4})


def test_tag_and_annotate():
    manifest = _many(3)
    manifest.tag("all", add=["raw"])
    manifest.tag(["c0", "c1"], add=["fine", "raw"], remove=())
    assert manifest["c0"].tags == ("raw", "fine")
    manifest.tag("c0", remove="raw")
    assert manifest["c0"].tags == ("fine",)
    manifest.annotate("c1", notes="odd case", metadata={"Re": 250})
    manifest.annotate("c1", metadata={"Ma": 0.1}, drop_metadata="Re")
    entry = manifest["c1"]
    assert entry.notes == "odd case" and entry.metadata == {"Ma": 0.1}
    assert entry.tags == ("raw", "fine")  # untouched by annotate


# --------------------------------------------------------------------------- #
# query filters                                                               #
# --------------------------------------------------------------------------- #
def test_entries_filters():
    manifest = DatasetManifest()
    manifest.add(
        "a.vtu",
        id="a",
        split="train",
        tags=["x", "y"],
        group="g/h",
        validate_source=False,
    )
    manifest.add(
        "b.vtu", id="b", split="test", tags=["x"], group="g", validate_source=False
    )
    manifest.add("c.vtu", id="c", validate_source=False)
    assert [e.id for e in manifest.entries()] == ["a", "b", "c"]
    assert [e.id for e in manifest.entries(split="train")] == ["a"]
    # an entry with no split matches only the no-filter call
    assert [e.id for e in manifest.entries(split="test")] == ["b"]
    assert [e.id for e in manifest.entries(tags=["x"])] == ["a", "b"]
    assert [e.id for e in manifest.entries(tags=["x", "y"])] == ["a"]
    assert [e.id for e in manifest.entries(group="g")] == ["a", "b"]
    assert [e.id for e in manifest.entries(group="g/h")] == ["a"]
    # path segments, not string prefixes
    manifest.annotate("c", group="gh")
    assert [e.id for e in manifest.entries(group="g")] == ["a", "b"]


def test_hand_edit_and_api_edit_interleave(tmp_path):
    """The single-source-of-truth rule: a hand edit between two API edits
    survives, because every API edit re-reads nothing and writes the whole
    document it holds -- so the workflow is load -> mutate -> save."""
    path = tmp_path / "m.json"
    manifest = DatasetManifest()
    manifest.add("a.vtu", id="a", validate_source=False)
    manifest.save(path)
    # hand edit: add a note
    doc = json.loads(path.read_text(encoding="utf-8"))
    doc["Entries"][0]["Notes"] = "hand-written"
    path.write_text(json.dumps(doc, indent=2) + "\n", encoding="utf-8")
    # tool edit: reload, curate, save
    manifest = DatasetManifest.load(path)
    manifest.set_split("a", "train")
    manifest.save(path)
    final = DatasetManifest.load(path)
    assert final["a"].notes == "hand-written"
    assert final["a"].split == "train"


def test_public_api_is_exported():
    for name in ("DatasetManifest", "DatasetEntry"):
        assert name in meshioplusplus.__all__
        assert hasattr(meshioplusplus, name)
    assert isinstance(DatasetEntry(id="x", source={"Path": "p"}), DatasetEntry)


# --------------------------------------------------------------------------- #
# paired coarse/fine entries                                                  #
# --------------------------------------------------------------------------- #
def _series(tmp_path, name, times, n=2):
    """A directory of single-step meshes whose filenames carry the step index."""
    import meshioplusplus

    directory = tmp_path / name
    directory.mkdir()
    for index, _ in enumerate(times):
        meshioplusplus.write(
            directory / f"out_{index:04d}.vtu", meshioplusplus.grid((n, n, n))
        )
    return directory


def test_target_is_optional_and_absent_serializes_identically(tmp_path):
    """The whole backward-compatibility claim: an entry with no Target must
    produce exactly the document it did before the key existed."""
    manifest = DatasetManifest()
    manifest.add({"Path": "a.vtu"}, validate_source=False)
    entry_doc = manifest.to_dict()["Entries"][0]
    assert set(entry_doc) == {"Id", "Source"}
    assert manifest["a"].target is None


def test_a_target_survives_the_save_round_trip(tmp_path):
    """`_with` rebuilds every entry on save; a field it forgets is dropped
    silently, so this is the guard for that whole class of bug."""
    path = tmp_path / "m.json"
    manifest = DatasetManifest()
    manifest.add(
        {"Pattern": "coarse/*.vtu"},
        target={"Pattern": "fine/*.vtu", "Format": "vtu"},
        id="sr",
        validate_source=False,
    )
    manifest.save(path)
    reloaded = DatasetManifest.load(path)
    assert reloaded["sr"].target == {"Pattern": "fine/*.vtu", "Format": "vtu"}
    # and a second save is byte-stable
    before = path.read_text(encoding="utf-8")
    reloaded.save(path)
    assert path.read_text(encoding="utf-8") == before
    assert "Target" in before


def test_target_time_series_and_entries_mirror_the_source_pair(tmp_path):
    coarse = _series(tmp_path, "coarse", [0, 1, 2], n=2)
    fine = _series(tmp_path, "fine", [0, 1, 2], n=4)
    manifest = DatasetManifest(base_dir=str(tmp_path))
    entry = manifest.add(
        {"Pattern": "coarse/*.vtu"}, target={"Pattern": "fine/*.vtu"}, id="sr"
    )
    assert len(entry.entries()) == len(entry.target_entries()) == 3
    assert len(entry.time_series()) == len(entry.target_time_series()) == 3
    _, coarse_mesh = entry.time_series()[0]
    _, fine_mesh = entry.target_time_series()[0]
    assert len(fine_mesh.points) > len(coarse_mesh.points)
    assert str(coarse) and str(fine)  # both directories were used


def test_an_entry_without_a_target_says_so_by_name(tmp_path):
    manifest = DatasetManifest()
    entry = manifest.add({"Path": "a.vtu"}, validate_source=False)
    with pytest.raises(ValueError, match="has no Target"):
        entry.target_time_series()
    with pytest.raises(ValueError, match="self-supervised"):
        entry.target_entries()


def test_pairing_refuses_a_mismatched_step_count(tmp_path):
    _series(tmp_path, "coarse", [0, 1, 2])
    _series(tmp_path, "fine", [0, 1])
    manifest = DatasetManifest(base_dir=str(tmp_path))
    with pytest.raises(ValueError, match="pairs 3 source steps with 2 target steps"):
        manifest.add(
            {"Pattern": "coarse/*.vtu"}, target={"Pattern": "fine/*.vtu"}, id="sr"
        )


def test_pairing_refuses_runs_at_different_instants(tmp_path):
    """The silent corruption the count check cannot see: two plans of the same
    length whose steps are at different times. Index pairing would train the
    model to map one moment onto another."""
    _series(tmp_path, "coarse", [0, 1, 2])
    _series(tmp_path, "fine", [0, 1, 2])
    manifest = DatasetManifest(base_dir=str(tmp_path))
    with pytest.raises(ValueError, match="are not at the same instants"):
        manifest.add(
            {"Pattern": "coarse/*.vtu", "Times": [0.0, 0.5, 1.0]},
            target={"Pattern": "fine/*.vtu", "Times": [0.0, 0.25, 0.5]},
            id="sr",
        )


def test_pairing_accepts_times_declared_meaningless(tmp_path):
    """`TimeFrom: index` says the times carry no information, so there is
    nothing for the two sides to disagree about -- the escape hatch the error
    message names.

    An explicit ``Times`` list always overrides ``TimeFrom`` (the plan then
    reports ``time_source="explicit"``), so the two cannot be combined to mean
    "ignore these times" -- declaring times is taken to mean them.
    """
    _series(tmp_path, "coarse", [0, 1, 2])
    _series(tmp_path, "fine", [0, 1, 2])
    manifest = DatasetManifest(base_dir=str(tmp_path))
    entry = manifest.add(
        {"Pattern": "coarse/*.vtu", "TimeFrom": "index"},
        target={"Pattern": "fine/*.vtu", "TimeFrom": "index"},
        id="sr",
    )
    assert check_pairing(entry) == 3

    with pytest.raises(ValueError, match="are not at the same instants"):
        DatasetManifest(base_dir=str(tmp_path)).add(
            {"Pattern": "coarse/*.vtu", "Times": [0.0, 0.5, 1.0], "TimeFrom": "index"},
            target={
                "Pattern": "fine/*.vtu",
                "Times": [0.0, 0.25, 0.5],
                "TimeFrom": "index",
            },
            id="sr",
        )


def test_check_pairing_counts_a_self_supervised_entry(tmp_path):
    _series(tmp_path, "coarse", [0, 1, 2])
    manifest = DatasetManifest(base_dir=str(tmp_path))
    entry = manifest.add({"Pattern": "coarse/*.vtu"}, id="solo")
    assert check_pairing(entry) == 3
