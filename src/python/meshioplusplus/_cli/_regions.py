import json

from .._helpers import read_metadata, reader_map


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
    parser.add_argument("--json", action="store_true", help="emit the regions as JSON")


def regions_cmd(args):
    # Goes through read_metadata rather than a full read: a native metadata
    # path costs nothing extra to report regions from an already-read mesh
    # (Exodus, and every fallback path), and this stays cheap on any format --
    # nothing here needs the connectivity, only the region list.
    meta = read_metadata(args.infile, file_format=args.input_format)
    regions = meta.get("regions") or []

    if args.json:
        print(json.dumps(regions, indent=2))
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
