"""tools/boolean_corpus.py (the surface-boolean spike, doc/spikes.md): the
case classes, the seeded pairs, the placement of a pair and the summary.
No network and no harness: the rows are written here."""

import importlib.util
import json
import pathlib

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[2]
_spec = importlib.util.spec_from_file_location(
    "boolean_corpus", ROOT / "tools" / "boolean_corpus.py"
)
corpus = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(corpus)


def test_a_case_name_gives_its_class():
    assert corpus.classify("seed78_12_sphere_box") == "sphere_box"
    assert corpus.classify("seed1_0_tori_crossing") == "tori_crossing"
    assert corpus.classify("thingi_self_123") == "thingi_self"
    assert corpus.classify("thingi_cross_1_2") == "thingi_cross"
    assert corpus.classify("boxes_shared_face") == "boxes_shared_face"


def test_the_pairs_are_seeded_and_cover_self_and_cross():
    a = corpus.pairs([3, 5, 9], 7808)
    assert a == corpus.pairs([3, 5, 9], 7808)
    assert a != corpus.pairs([3, 5, 9], 1)
    assert [p["case"] for p in a] == [
        "thingi_self_3",
        "thingi_cross_3_5",
        "thingi_self_5",
        "thingi_cross_5_9",
        "thingi_self_9",
        "thingi_cross_9_3",
    ]
    for p in a:
        r = np.asarray(p["rotation"])
        assert np.allclose(r @ r.T, np.eye(3)) and np.isclose(np.linalg.det(r), 1.0)
        n = np.linalg.norm(p["shift"])
        assert (0.2 - 1e-12 <= n <= 0.5 + 1e-12) if p["b"] == p["a"] else n <= 0.33
    assert [p["case"] for p in corpus.pairs([4], 0)] == ["thingi_self_4"]


def test_placing_keeps_the_shape_and_moves_it_by_its_size():
    pts = np.array([[0, 0, 0], [2, 0, 0], [0, 4, 0], [0, 0, 1]], dtype=float)
    ident = np.eye(3).tolist()
    assert np.allclose(corpus.place(pts, ident, [0, 0, 0], False), pts)
    moved = corpus.place(pts, ident, [0.5, 0, 0], False)
    assert np.allclose(moved - pts, [[2.0, 0, 0]] * 4)  # half of the size, 4
    unit = corpus.place(pts, ident, [0, 0, 0], True)
    assert np.isclose(np.ptp(unit, axis=0).max(), 1.0)
    assert np.allclose(unit.min(axis=0) + unit.max(axis=0), 0.0)
    # A rotation keeps every distance.
    rot = corpus.pairs([1], 3)[0]["rotation"]
    turned = corpus.place(pts, rot, [0, 0, 0], False)
    d = np.linalg.norm(pts[:, None] - pts[None], axis=2)
    assert np.allclose(np.linalg.norm(turned[:, None] - turned[None], axis=2), d)


def test_the_summary_counts_each_class(tmp_path):
    rows = [
        {
            "case": f"seed1_{i}_spheres_offset",
            "op": "union",
            "ok": True,
            "watertight": True,
            "self_intersection_free": i != 0,
            "rounding_inverted": 1 if i == 0 else 0,
            "identity_error_union": 1e-15,
            "identity_error_difference": 2e-15,
            "seconds": 0.1 * (i + 1),
        }
        for i in range(2)
    ]
    rows.append(
        json.loads(corpus._failed("boxes_shared_face", "union", "coplanar contact"))
    )
    table = corpus.render(rows)
    assert (
        "| spheres_offset | 2 | 100% | 100% | 50% | 50% | 2.0e-15 | 0.150 / 0.200 | - |"
        in table
    )
    assert "| boxes_shared_face | 1 | 0% |" in table and "coplanar contact (1)" in table
    assert table.endswith(
        "3 runs: 2 succeeded, 1 watertight and free of self-intersections."
    )
    runs = tmp_path / "runs.jsonl"
    runs.write_text("\n".join(json.dumps(r) for r in rows) + "\n")
    assert corpus.main(["summary", str(runs)]) == 0
