import numpy as np

from .._common import warn
from .._helpers import read, read_metadata, reader_map


def add_args(parser):
    parser.add_argument("infile", type=str, help="mesh file to be read from")
    parser.add_argument(
        "--input-format",
        "-i",
        type=str,
        choices=sorted(list(reader_map.keys())),
        help="input file format",
        default=None,
    )
    parser.add_argument(
        "--fast",
        action="store_true",
        help=(
            "summarize from the file header instead of loading it. Skips "
            "decoding the data arrays, so it is much cheaper on large files "
            "(dramatically so for XDMF, whose DataItems declare their shape); "
            "formats with no header-only path are still read in full and say "
            "so. Omits the point/cell consistency checks, which need the "
            "connectivity."
        ),
    )


def _print_metadata(meta):
    print("<meshio++ mesh summary>")
    print(f"  Format: {meta['format'] or 'unknown'}")
    print(f"  Number of points: {meta['num_points']}")
    if meta["cell_blocks"]:
        print("  Number of cells:")
        for block in meta["cell_blocks"]:
            print(f"    {block['type']}: {block['num_cells']}")
    else:
        print("  No cells.")

    for label, key in (
        ("Point data", "point_data_names"),
        ("Cell data", "cell_data_names"),
        ("Field data", "field_data_names"),
    ):
        if meta[key]:
            print(f"  {label}: {', '.join(meta[key])}")

    # Only worth printing when there is a choice to make: a single-step file
    # gives a caller nothing to pass to --time-step.
    times = meta.get("time_values") or []
    if len(times) > 1:
        shown = ", ".join(f"{v:g}" for v in times[:8])
        if len(times) > 8:
            shown += ", ..."
        print(f"  Time steps: {len(times)} [{shown}]")

    regions = meta.get("regions") or []
    if regions:
        print(f"  Regions ({len(regions)}):")
        for r in regions:
            tag = "" if r["tag"] < 0 else f" tag={r['tag']}"
            print(f"    {r['name']} ({r['kind']}, {r['num_entries']} entries{tag})")

    provenance = meta.get("provenance") or []
    if provenance:
        # Say whose block it is: a file can carry a comment that merely looks
        # like a header, and reporting that as meshio++'s would be a lie.
        whose = (
            "" if meta.get("provenance_recognised") else " (not written by meshio++)"
        )
        print(f"  Provenance{whose}:")
        for line in provenance:
            print(f"    {line}")

    if "bbox_min" in meta:
        lo = ", ".join(f"{v:g}" for v in meta["bbox_min"])
        hi = ", ".join(f"{v:g}" for v in meta["bbox_max"])
        print(f"  Bounding box: [{lo}] - [{hi}]")

    # Say plainly when "fast" wasn't, rather than implying a saving that did
    # not happen.
    if meta["fell_back_to_full_read"]:
        print("  (no header-only path for this format; the file was read in full)")


def info(args):
    if args.fast:
        _print_metadata(read_metadata(args.infile, file_format=args.input_format))
        return 0

    # read mesh data
    mesh = read(args.infile, file_format=args.input_format)
    print(mesh)

    # check if the cell arrays are consistent with the points
    is_consistent = True
    for cells in mesh.cells:
        if np.any(cells.data > mesh.points.shape[0]):
            warn("Inconsistent mesh. Cells refer to nonexistent points.")
            is_consistent = False
            break

    # check if there are redundant points
    if is_consistent:
        point_is_used = np.zeros(mesh.points.shape[0], dtype=bool)
        for cells in mesh.cells:
            point_is_used[cells.data] = True
        if np.any(~point_is_used):
            warn("Some points are not part of any cell.")

    return 0
