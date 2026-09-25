# MIT License -- part of meshio++ (https://github.com/loumalouomega/meshioplusplus).
"""Export a Femap model (``.modfem``) as a neutral file (``.neu``) meshio++ reads.

``.modfem`` is Femap's private database; the only way out is Femap itself. This
script drives a running (or newly started) Femap on Windows through its COM API
with ``pywin32`` and asks it to write the model as a neutral file, which
meshio++ reads natively (``doc/formats/femap.md``)::

    python export_neutral.py model.modfem [model.neu] [--version 2401]

It opens the model, writes nodes, elements, properties, materials, groups and
output sets, and leaves Femap running. The API call is ``feFileWriteNeutral``
(Femap API reference, "File" methods); its argument list is the one to check
against your Femap version's API help first (see ``doc/routes/femap_modfem.md``).
"""

from __future__ import print_function

import argparse
import os

FE_OK = -1  # Femap's success code


def export(femap, model, out, version):
    rc = femap.feFileOpen(True, os.path.abspath(model))
    if rc != FE_OK:
        raise SystemExit(
            "export_neutral: Femap could not open %s (code %s)" % (model, rc)
        )
    # Group 0 (the whole model); geometry off; model, groups and output on.
    rc = femap.feFileWriteNeutral(
        0, os.path.abspath(out), False, True, True, True, False, float(version), 0
    )
    if rc != FE_OK:
        raise SystemExit(
            "export_neutral: Femap could not write %s (code %s)" % (out, rc)
        )


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("model")
    parser.add_argument(
        "out", nargs="?", help="the .neu to write (default: <model>.neu)"
    )
    parser.add_argument(
        "--version",
        default="2401",
        help="neutral file version to write (default: %(default)s)",
    )
    args = parser.parse_args(argv)
    out = args.out or os.path.splitext(args.model)[0] + ".neu"

    import win32com.client

    try:
        femap = win32com.client.GetActiveObject("femap.model")
    except Exception:  # no running Femap: start one
        femap = win32com.client.Dispatch("femap.model")
    export(femap, args.model, out, args.version)
    print("export_neutral: wrote %s" % out)


if __name__ == "__main__":
    main()
