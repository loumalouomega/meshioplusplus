import json

from .._helpers import read, reader_map
from .._stats import compute_stats


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
        "--json", action="store_true", help="emit the statistics as JSON"
    )


def _fmt3(v):
    return "[" + ", ".join(f"{x:.6g}" for x in v) + "]"


def stats_cmd(args):
    mesh = read(args.infile, file_format=args.input_format)
    s = compute_stats(mesh)

    if args.json:
        print(json.dumps(s, indent=2))
        return 0

    print("<meshio++ geometric stats>")
    print(f"  points:          {s['num_points']}")
    print(f"  cells:           {s['num_cells']}")
    print(f"  bbox min:        {_fmt3(s['bbox_min'])}")
    print(f"  bbox max:        {_fmt3(s['bbox_max'])}")
    print(f"  extent:          {_fmt3(s['extent'])}")
    print(f"  centroid:        {_fmt3(s['centroid'])}")
    print("  cell types:")
    for name, count in s["cell_type_counts"].items():
        print(f"    {name}: {count}")
    print(f"  total area:      {s['total_area']:.6g}")
    print(f"  signed volume:   {s['signed_volume']:.6g}")
    print(f"  unsigned volume: {s['unsigned_volume']:.6g}")
    print(f"  inverted cells:  {s['num_inverted']}")
    return 0
