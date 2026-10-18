"""Freeze the column sets of MSC Nastran's HDF5 tables from the committed fixtures.

The writers of ``nastran_h5`` give every table they write MSC's own columns, in MSC's
order, with MSC's types and array lengths, so that a reader of the schema finds what it
expects. This reads the six files of ``tests/python/meshes/nastran_h5`` and writes the
column sets of the ``INPUT`` element, property and material tables and the ``RESULT``
nodal and elemental ones, as text, to

- ``src/python/meshioplusplus/nastran_h5/_schema.py`` (a dict), and
- ``src/cpp/src/formats/nastran_h5.cpp`` (the block between the ``BEGIN/END GENERATED
  SCHEMA`` markers),

both parsed at run time by the writers the same way. A column reads ``NAME:i8``, ``f8``,
``S4`` or ``S8``, with ``[n]`` for an array of ``n``.

    python tools/gen_nastran_h5_schema.py
"""

import glob
import pathlib

import h5py

ROOT = pathlib.Path(__file__).resolve().parent.parent
FIXTURES = sorted(glob.glob(str(ROOT / "tests/python/meshes/nastran_h5/*.h5")))
ROOTS = (
    "INPUT/ELEMENT",
    "INPUT/PROPERTY",
    "INPUT/MATERIAL",
    "RESULT/NODAL",
    "RESULT/ELEMENTAL",
)
KINDS = {"int64": "i8", "float64": "f8", "|S4": "S4", "|S8": "S8"}


def column(dtype, name):
    base = dtype[name].base
    kind = KINDS[str(base)]
    shape = dtype[name].shape
    return f"{name}:{kind}" + (f"[{shape[0]}]" if shape else "")


def tables(f):
    found = {}

    def visit(name, item):
        if isinstance(item, h5py.Dataset) and item.dtype.names:
            if any(name.startswith("NASTRAN/" + r + "/") for r in ROOTS):
                found[name[len("NASTRAN/") :]] = " ".join(
                    column(item.dtype, n) for n in item.dtype.names
                )

    f.visititems(visit)
    return found


def main():
    schema = {}
    for fixture in FIXTURES:
        with h5py.File(fixture) as f:
            for path, text in tables(f).items():
                if schema.setdefault(path, text) != text:
                    raise SystemExit(
                        f"{path} differs between fixtures:\n{schema[path]}\n{text}"
                    )
    lines = [f'    "{path}": "{schema[path]}",' for path in sorted(schema)]
    py = (
        '"""MSC Nastran HDF5 column sets, frozen from the fixtures by\n'
        "``tools/gen_nastran_h5_schema.py``. Do not edit.\n\n"
        "Each table's columns as text: ``NAME:i8``, ``f8``, ``S4`` or ``S8``, with ``[n]`` for\n"
        'an array of ``n``."""\n\n'
        "SCHEMA = {\n" + "\n".join(lines) + "\n}\n"
    )
    (ROOT / "src/python/meshioplusplus/nastran_h5/_schema.py").write_text(py)
    cpp_lines = [f'    {{"{path}", "{schema[path]}"}},' for path in sorted(schema)]
    block = (
        "// BEGIN GENERATED SCHEMA (tools/gen_nastran_h5_schema.py)\n"
        "constexpr Nh5wSchemaEntry kNh5wSchema[] = {\n"
        + "\n".join(cpp_lines)
        + "\n};\n"
        "// END GENERATED SCHEMA\n"
    )
    path = ROOT / "src/cpp/src/formats/nastran_h5.cpp"
    text = path.read_text()
    begin = text.find("// BEGIN GENERATED SCHEMA")
    end = text.find("// END GENERATED SCHEMA\n")
    if begin < 0 or end < 0:
        raise SystemExit("the markers are missing from nastran_h5.cpp")
    path.write_text(
        text[:begin] + block + text[end + len("// END GENERATED SCHEMA\n") :]
    )
    print(f"{len(schema)} tables")


if __name__ == "__main__":
    main()
