"""Unit tests for the CLI's shared argument helpers (``_cli/_common.py``).

The verbs' ``--help`` output is pinned by ``test_cli.py``; these tests pin the
helpers themselves: the value parsers and the option blocks they add.
"""

import argparse

import pytest

from meshioplusplus._cli import _common


def test_split_csv_strips_and_drops_empty_entries():
    assert _common.split_csv("a, b,,c ,") == ["a", "b", "c"]
    assert _common.split_csv("") == []
    assert _common.split_csv(None) == []
    assert _common.split_csv(" , ") == []


@pytest.mark.parametrize("text", [None])
def test_parse_helpers_pass_none_through(text):
    assert _common.parse_bounds(text, "verb") is None
    assert _common.parse_ints(text, "verb", "--resolution") is None


def test_parse_bounds():
    assert _common.parse_bounds("0,1,2,3,4,5", "sdf") == [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
    assert _common.parse_bounds("-1,-1,-1,1,1,1", "sdf")[0] == -1.0
    with pytest.raises(ValueError, match=r"^sdf: --bounds expects 'xlo,ylo,zlo"):
        _common.parse_bounds("0,1,2", "sdf")
    with pytest.raises(ValueError):
        _common.parse_bounds("0,1,2,3,4,x", "sdf")


def test_parse_ints_counts_and_message():
    assert _common.parse_ints("4,5,6", "sdf", "--resolution") == [4, 5, 6]
    with pytest.raises(ValueError, match=r"^sdf: --resolution expects 'nx,ny,nz'$"):
        _common.parse_ints("4,5", "sdf", "--resolution")
    assert _common.parse_ints("4,5", "voxelize", "--resolution", counts=(2, 3)) == [
        4,
        5,
    ]
    with pytest.raises(ValueError, match=r"\(or 'nx,ny'\)$"):
        _common.parse_ints(
            "4",
            "voxelize",
            "--resolution",
            counts=(2, 3),
            form="'nx,ny,nz' (or 'nx,ny')",
        )
    with pytest.raises(ValueError):
        _common.parse_ints("4,a,6", "sdf", "--resolution")


def _dests(parser):
    return [a.dest for a in parser._actions if a.dest != "help"]


def test_add_box_args_order_defaults_and_help():
    parser = argparse.ArgumentParser()
    _common.add_box_args(parser, bounds_of="skin's", padding_relative=0.1)
    assert _dests(parser) == ["bounds", "padding", "padding_relative"]
    args = parser.parse_args([])
    assert (args.bounds, args.padding, args.padding_relative) == (None, 0.0, 0.1)
    text = " ".join(parser.format_help().split())
    assert "the skin's bounding box by default" in text
    custom = argparse.ArgumentParser()
    _common.add_box_args(custom, bounds_help="only this box")
    assert "only this box" in custom.format_help()
    assert custom.parse_args([]).padding_relative == 0.0


def test_add_transfer_io_args_and_flags():
    parser = argparse.ArgumentParser()
    _common.add_transfer_io_args(parser, "out help")
    _common.add_on_conflict_arg(parser)
    _common.add_transfer_quiet_arg(parser)
    assert _dests(parser) == [
        "source",
        "target",
        "outfile",
        "input_format",
        "output_format",
        "on_conflict",
        "quiet",
    ]
    args = parser.parse_args(["a.vtu", "b.vtu", "c.vtu", "-q"])
    assert args.on_conflict == "error" and args.quiet
    with pytest.raises(SystemExit):
        parser.parse_args(["a", "b", "c", "--on-conflict", "nope"])


def test_read_source_target(tmp_path, monkeypatch):
    seen = []

    def fake_read(path, file_format=None):
        seen.append((path, file_format))
        return path

    monkeypatch.setattr("meshioplusplus._helpers.read", fake_read)
    args = argparse.Namespace(source="s", target="t", input_format="vtu")
    assert _common.read_source_target(args) == ("s", "t")
    assert seen == [("s", "vtu"), ("t", "vtu")]


_HELP = {
    "in": "in",
    "ratio": "ratio",
    "count": "count",
    "max_error": "err",
    "placement": "place",
    "features": "feat",
    "angle": "angle",
}


def test_add_decimate_args_block():
    parser = argparse.ArgumentParser()
    boundary = (("--no-preserve-boundary",), {"dest": "b", "action": "store_false"})
    _common.add_decimate_args(parser, "faces", _HELP, boundary)
    assert _dests(parser) == [
        "infile",
        "outfile",
        "input_format",
        "output_format",
        "ratio",
        "target_faces",
        "max_error",
        "placement",
        "b",
        "preserve_features",
        "feature_angle",
        "quiet",
    ]
    args = parser.parse_args(["a", "b", "--target-faces", "10"])
    assert args.target_faces == 10 and args.b and args.feature_angle == 30.0
    # The stop criteria are one required, mutually exclusive group.
    with pytest.raises(SystemExit):
        parser.parse_args(["a", "b"])
    with pytest.raises(SystemExit):
        parser.parse_args(["a", "b", "--ratio", "0.5", "--max-error", "1"])


def test_decimate_summary():
    class Block:
        data = [0, 1, 2]

    class Out:
        cells = [Block(), Block()]

    report = {
        "tets_removed": 4,
        "points_removed": 2,
        "collapses_rejected": 1,
        "max_error_applied": 0.125,
    }
    assert _common.decimate_summary(Out(), report, "tets", "tets_removed") == (
        "decimated to 6 tets (4 removed, 2 points removed, "
        "1 collapses rejected, max error 0.125)"
    )
