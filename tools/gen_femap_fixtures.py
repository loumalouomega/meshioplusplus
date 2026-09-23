#!/usr/bin/env python3
"""Regenerate the Femap neutral fixtures under ``tests/python/meshes/femap/``.

Every mesh here was written by Femap or by a program that exports to it, never
by meshio++, and all of it is MIT-licensed (see that directory's README.md):

    python tools/gen_femap_fixtures.py FRONTISTR_EXAMPLES EMS_REPO NEUTRAL_PARSER_REPO

* ``FRONTISTR_EXAMPLES`` is ``fistr1/tools/neu2fstr/example`` of FrontISTR
  (https://github.com/FrontISTR/FrontISTR): real Femap 8.2 models. Copied
  unmodified: ``A342`` (tetra10), ``A352`` (wedge15), ``A362`` (brick20),
  ``A731`` (tri3), ``B741`` (quad4) and ``MC361`` (brick8 with node groups).
* ``EMS_REPO`` is https://github.com/EMSolution-SSIL/ems_file_format_converter:
  Femap 4.41 files EMSolution writes. ``mesh_sample.neu`` (tetra, pyramid, wedge,
  brick) is copied; ``post_sample451.neu`` and ``post_sample1051.neu`` hold 13
  output sets on that same mesh, so they are appended to it as
  ``ems_results_451.neu`` and ``ems_results_1051.neu`` -- one file per results
  encoding, holding the same values.
* ``NEUTRAL_PARSER_REPO`` is https://framagit.org/numenic/femap_neutral_parser:
  results-only files written by Femap 2020.1, Femap 8.2 and MYSTRAN for a
  12-node, 11-bar model whose geometry they do not carry. Each is prefixed with
  that model (``403``/``404`` blocks written here: nodes 1-12 along x, bar ``i``
  joining nodes ``i`` and ``i+1``) as ``v2020_results.neu``, ``v82_results.neu``
  and ``mystran_results.neu``. The header block of the original is kept, so the
  version the result blocks were written in is the one the file declares.
"""

import pathlib
import shutil
import sys

OUT = (
    pathlib.Path(__file__).resolve().parent.parent
    / "tests"
    / "python"
    / "meshes"
    / "femap"
)


def read_text(path):
    return path.read_bytes().decode("latin-1").replace("\r\n", "\n")


def split_header(text):
    """(everything up to and including the 100 block, the rest)."""
    lines = text.split("\n")
    for i, line in enumerate(lines):
        if (
            line.strip() == "-1"
            and i + 1 < len(lines)
            and lines[i + 1].strip() == "100"
        ):
            j = i + 2
            while lines[j].strip() != "-1":
                j += 1
            return "\n".join(lines[: j + 1]) + "\n", "\n".join(lines[j + 1 :])
    raise SystemExit("no header block")


def bar_model():
    nodes = "".join(
        f"{i},0,0,1,46,0,0,0,0,0,0,{float(i - 1)},0.,0.,0,\n" for i in range(1, 13)
    )
    elements = []
    for i in range(1, 12):
        slots = [i, i + 1] + [0] * 18
        elements.append(
            f"{i},124,1,1,0,1,0,0,0,0,0,0,0,\n"
            + "".join(f"{v}," for v in slots[:10])
            + "\n"
            + "".join(f"{v}," for v in slots[10:])
            + "\n0.,0.,0.,\n0.,0.,0.,\n0.,0.,0.,\n0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,\n"
        )
    return (
        "   -1\n   403\n"
        + nodes
        + "   -1\n   -1\n   404\n"
        + "".join(elements)
        + "   -1\n"
    )


def main():
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    frontistr, ems, parser = (pathlib.Path(p) for p in sys.argv[1:])
    OUT.mkdir(parents=True, exist_ok=True)
    for name in ["A/A342", "A/A352", "A/A362", "A/A731", "B/B741", "heat/MC361"]:
        src = frontistr / f"{name}.NEU"
        shutil.copyfile(src, OUT / f"{src.stem}.neu")
    data = ems / "sample"
    if not (data / "mesh_sample.neu").exists():
        data = ems
    shutil.copyfile(data / "mesh_sample.neu", OUT / "mesh_sample.neu")
    mesh = read_text(data / "mesh_sample.neu")
    for kind in ("451", "1051"):
        _, results = split_header(read_text(data / f"post_sample{kind}.neu"))
        (OUT / f"ems_results_{kind}.neu").write_text(
            mesh.rstrip("\n") + "\n" + results, newline="\n"
        )
    samples = parser / "tests" / "data"
    if not samples.exists():
        samples = parser
    for src, dst in [
        ("FEMAP_v2020-1-0.neu", "v2020_results.neu"),
        ("FEMAP_v8-2.neu", "v82_results.neu"),
        ("mystran_00.NEU", "mystran_results.neu"),
    ]:
        header, results = split_header(read_text(samples / src))
        (OUT / dst).write_text(header + bar_model() + results, newline="\n")


if __name__ == "__main__":
    main()
