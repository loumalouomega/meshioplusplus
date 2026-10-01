import csv
import sys

from .._helpers import read, reader_map
from .._interfaces import contact_pairs
from ._json import emit_json


def add_args(parser):
    parser.add_argument(
        "slave_file", type=str, help="mesh containing the slave Point region"
    )
    parser.add_argument(
        "master_file", type=str, help="mesh containing the master Cell region"
    )
    parser.add_argument("--slave-region", required=True, help="Point region to project")
    parser.add_argument("--master-region", required=True, help="master Cell region")
    parser.add_argument("--slave-format", choices=sorted(reader_map), default=None)
    parser.add_argument("--master-format", choices=sorted(reader_map), default=None)
    parser.add_argument(
        "--output", "-o", default=None, help="optional CSV file for the pair table"
    )
    parser.add_argument("--tolerance", type=float, default=0.0)
    parser.add_argument("--require-complete", action="store_true")
    parser.add_argument(
        "--json", action="store_true", help="emit the full pair table as JSON"
    )


def contact_pairs_cmd(args):
    slave = read(args.slave_file, file_format=args.slave_format)
    master = read(args.master_file, file_format=args.master_format)
    result = contact_pairs(
        slave,
        args.slave_region,
        args.master_region,
        master_mesh=master,
        tolerance=args.tolerance,
        require_complete=args.require_complete,
    )
    pairs = []
    for i in range(len(result["slave_point"])):
        pairs.append(
            {
                "slave_point": int(result["slave_point"][i]),
                "master_cell": int(result["master_cell"][i]),
                "master_facet": int(result["master_facet"][i]),
                "master_subfacet": int(result["master_subfacet"][i]),
                "local_coordinates": result["local_coordinates"][i].tolist(),
                "closest_point": result["closest_point"][i].tolist(),
                "gap": float(result["gap"][i]),
                "normal": result["normal"][i].tolist(),
            }
        )
    if args.json:
        emit_json(
            {
                "num_pairs": len(pairs),
                "pairs": pairs,
                "unmatched": result["unmatched"].tolist(),
            }
        )
    elif args.output:
        with open(args.output, "w", newline="", encoding="utf-8") as stream:
            writer = csv.writer(stream)
            writer.writerow(
                [
                    "slave_point",
                    "master_cell",
                    "master_facet",
                    "master_subfacet",
                    "u",
                    "v",
                    "w",
                    "x",
                    "y",
                    "z",
                    "gap",
                    "nx",
                    "ny",
                    "nz",
                ]
            )
            for pair in pairs:
                writer.writerow(
                    [
                        pair["slave_point"],
                        pair["master_cell"],
                        pair["master_facet"],
                        pair["master_subfacet"],
                        *pair["local_coordinates"],
                        *pair["closest_point"],
                        pair["gap"],
                        *pair["normal"],
                    ]
                )
        print(
            f"contact pairs: {len(pairs)} node(s), {len(result['unmatched'])} unmatched; wrote {args.output}"
        )
    else:
        writer = csv.writer(sys.stdout)
        writer.writerow(
            ["slave_point", "master_cell", "master_facet", "master_subfacet", "gap"]
        )
        writer.writerows(
            [
                [
                    p["slave_point"],
                    p["master_cell"],
                    p["master_facet"],
                    p["master_subfacet"],
                    p["gap"],
                ]
                for p in pairs
            ]
        )
    return 0
