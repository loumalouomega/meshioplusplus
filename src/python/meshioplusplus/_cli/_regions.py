import argparse

from .._helpers import _writer_map, read, read_metadata, reader_map, write
from .._region_ops import edit_regions, parse_edit
from ._json import emit_json

_EDIT_FLAGS = (
    ("union", "OUT=A,B", "add OUT, the union of regions A, B, ..."),
    ("intersection", "OUT=A,B", "add OUT, the entries in every one of A, B, ..."),
    ("difference", "OUT=A,B", "add OUT, the entries of A in none of B, ..."),
    ("rename", "OLD=NEW", "rename a region"),
    ("retag", "NAME=TAG[:DIM]", "set a region's tag (and dimension)"),
    ("delete", "NAME", "remove a region"),
)


class _EditAction(argparse.Action):
    """Collect the edit flags in command-line order."""

    def __call__(self, parser, namespace, values, option_string=None):
        edits = list(getattr(namespace, "edits", None) or [])
        edits.append((self.dest.removeprefix("edit_"), values))
        namespace.edits = edits


def add_args(parser):
    parser.add_argument("infile", type=str, help="mesh file to be read from")
    parser.add_argument(
        "outfile",
        type=str,
        nargs="?",
        default=None,
        help="with edits: the edited mesh to be written to",
    )
    parser.add_argument(
        "--input-format",
        "-i",
        type=str,
        choices=sorted(list(reader_map.keys())),
        help="input file format",
        default=None,
    )
    parser.add_argument(
        "--output-format",
        "-o",
        type=str,
        choices=sorted(list(_writer_map.keys())),
        help="output file format",
        default=None,
    )
    parser.add_argument("--json", action="store_true", help="emit the regions as JSON")
    group = parser.add_argument_group(
        "edits",
        "applied in command-line order; prefix a name with point:, cell: or side: to "
        "pin its kind",
    )
    for op, metavar, what in _EDIT_FLAGS:
        group.add_argument(
            f"--{op}", dest=f"edit_{op}", metavar=metavar, action=_EditAction, help=what
        )
    group.add_argument(
        "--drop-inputs",
        action="store_true",
        help="set operations remove their input regions",
    )
    parser.set_defaults(edits=[])


def _describe(regions):
    return [
        {
            "name": r.name,
            "kind": r.kind,
            "dim": r.dim,
            "tag": r.tag,
            "num_entries": len(r.entries),
        }
        for r in regions
    ]


def _edit_cmd(args):
    if not args.outfile:
        raise SystemExit("regions: edits need an output file")
    mesh = read(args.infile, file_format=args.input_format)
    edits = [parse_edit(op, text, not args.drop_inputs) for op, text in args.edits]
    out = edit_regions(mesh, edits)
    write(args.outfile, out, file_format=args.output_format)
    regions = _describe(out.regions)
    if args.json:
        emit_json(regions)
        return 0
    print(f"<meshio++ mesh regions> ({len(regions)}) after {len(edits)} edit(s)")
    for r in regions:
        tag = "" if r["tag"] < 0 else f", tag={r['tag']}"
        dim = "" if r["dim"] < 0 else f", dim={r['dim']}"
        print(f"  {r['name']} ({r['kind']}, {r['num_entries']} entries{dim}{tag})")
    return 0


def regions_cmd(args):
    if args.edits:
        return _edit_cmd(args)
    if args.outfile:
        raise SystemExit("regions: an output file needs at least one edit")
    # Goes through read_metadata rather than a full read: a native metadata
    # path costs nothing extra to report regions from an already-read mesh
    # (Exodus, and every fallback path), and this stays cheap on any format --
    # nothing here needs the connectivity, only the region list.
    meta = read_metadata(args.infile, file_format=args.input_format)
    regions = meta.get("regions") or []

    if args.json:
        emit_json(regions)
        return 0

    if not regions:
        print("<meshio++ mesh regions>")
        print("  No regions.")
        if meta.get("fell_back_to_full_read"):
            print("  (the full mesh was read; this format may simply carry none)")
        return 0

    print(f"<meshio++ mesh regions> ({len(regions)})")
    for r in regions:
        tag = "" if r["tag"] < 0 else f", tag={r['tag']}"
        dim = "" if r["dim"] < 0 else f", dim={r['dim']}"
        print(f"  {r['name']} ({r['kind']}, {r['num_entries']} entries{dim}{tag})")
    return 0
