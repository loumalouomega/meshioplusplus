"""Write the export table of the Windows shared core (a .def file).

WINDOWS_EXPORT_ALL_SYMBOLS exports every external symbol of every object, the
standard library's, Eigen's and pugixml's instantiations among them, and that
is more than the linker's 65,535 exports (LNK1189). The add-on providers need
the core's own code only, so this keeps the symbols whose decorated name sits
in the ``meshioplusplus`` namespace and nothing else.

    gen_core_def.py --dumpbin DUMPBIN --out core.def OBJDIR [OBJDIR ...]

Run between compiling and linking (a PRE_LINK step), so it reads the objects
the build just wrote under each OBJDIR.
"""

import argparse
import pathlib
import subprocess
import sys

LIMIT = 65000  # the linker stops at 65,535 with room kept for the compiler's own


def external_symbols(dumpbin, obj):
    """Yield (name, is_data) for each external symbol the object defines."""
    out = subprocess.run(
        [dumpbin, "/NOLOGO", "/SYMBOLS", str(obj)],
        capture_output=True,
        text=True,
        check=True,
    ).stdout
    for line in out.splitlines():
        if " External " not in line or "|" not in line or "UNDEF" in line:
            continue
        head, _, tail = line.partition("|")
        if " SECT" not in head:
            continue
        name = tail.strip().split(" ", 1)[0]
        if not name or name[0] == "$":
            continue
        yield name, "()" not in head


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dumpbin", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--namespace", default="meshioplusplus@@")
    parser.add_argument("dirs", nargs="+")
    args = parser.parse_args()

    exports = {}
    count = 0
    for directory in args.dirs:
        for obj in sorted(pathlib.Path(directory).rglob("*.obj")):
            count += 1
            for name, is_data in external_symbols(args.dumpbin, obj):
                # String literals and the compiler's own helpers are never API.
                if args.namespace in name and not name.startswith("??_C"):
                    exports[name] = exports.get(name, False) or is_data
    if not count:
        sys.exit("gen_core_def: no .obj files under " + ", ".join(args.dirs))
    if len(exports) > LIMIT:
        sys.exit(f"gen_core_def: {len(exports)} exports, over the {LIMIT} budget")

    lines = ["EXPORTS"]
    lines += [f"    {n}{' DATA' if d else ''}" for n, d in sorted(exports.items())]
    pathlib.Path(args.out).write_text("\n".join(lines) + "\n")
    print(f"gen_core_def: {len(exports)} exports from {count} objects -> {args.out}")


if __name__ == "__main__":
    main()
