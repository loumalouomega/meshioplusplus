"""``detect INFILE [--json]``: the ranked format report (roadmap 2.4.9)."""

from .._sniff import detect
from ._json import emit_json


def add_args(parser):
    parser.add_argument("infile", type=str, help="file or directory to identify")
    parser.add_argument("--json", action="store_true", help="print the report as JSON")


def detect_cmd(args):
    found = detect(args.infile)
    if args.json:
        emit_json(
            {
                "path": args.infile,
                "candidates": [
                    {"format": f, "confidence": c, "reason": r} for f, c, r in found
                ],
            }
        )
    elif not found:
        print(f"no format matches {args.infile}")
    else:
        for fmt, confidence, reason in found:
            print(f"{fmt}  {confidence}  ({reason})")
    return 0 if found else 1
