# Documentation diagrams

Every figure under `doc/public/diagrams/` is **generated** by `gen_diagrams.py` — the architecture boxes-and-arrows, the flow diagrams and the cell-type node orderings alike — and committed together with a 2x PNG twin. Nothing in that directory is edited by hand.

```
doc/diagrams/
  gen_diagrams.py        the CLI: write everything, or --check / --list / --only NAME / --no-png
  diaglib/palette.py     the one place colours, fonts and stroke widths live
  diaglib/svg.py         a deterministic SVG builder (two-decimal coordinates, fixed attribute order)
  diaglib/tables.py      loaders for the code's own topology tables + reference-element corners
  diaglib/project.py     the orthographic camera (tikz-3dplot's 70/110 main coords) and depth sort
  diaglib/cellart.py     drawing a reference element: hidden edges dashed, numbered nodes
  figures/*.py           one function per figure, registered in figures/REGISTRY
```

## Regenerate and check

```sh
python3 doc/diagrams/gen_diagrams.py            # every SVG + PNG twin into doc/public/diagrams/
python3 doc/diagrams/gen_diagrams.py --check    # what CI's lint job runs
```

The generator is standard library only, so it runs on any Python 3.9+ with no build and no `meshioplusplus` import. PNG twins are rasterised at 2x with `rsvg-convert` (falling back to `inkscape`, then ImageMagick); their bytes are not pinned because rasterisers differ per machine, only that they exist at the expected size.

`--check` fails when a committed SVG differs from what the generator produces (stale), when a figure has no PNG twin, when a figure is embedded by no page under `doc/` nor by `README.md`, or when a page references a figure no function produces. `tests/python/test_doc_diagrams.py` pins the same things and, in addition, that every node of the cell-type figures sits where the code's tables say it does.

## Where the geometry comes from

The cell-type figures do not transcribe node orderings. `diaglib/tables.py` reads the literal tables straight from the package source with `ast` (`_skin._CELL_FACES`, `_surface._CELL_EDGES`, `_convert_cells._ELEVATE`, `_common.num_nodes_per_cell`) and loads `_refine_templates.py` as a module (it imports only `itertools`); mid-edge nodes are the midpoints those tables name and face centres are the centroids of `QUAD_FACES` rows. The only hand-typed geometry is each linear reference element's corner positions, taken from the TikZ file `doc/cell_types.tex` that these figures replaced (deleted in the same change; git history keeps it).

## Style rules

- Every figure sits on an explicit light paper rectangle. An `<img>`-loaded SVG cannot see VitePress's dark-mode class and a `prefers-color-scheme` query would follow the operating system rather than the site toggle, so figures are self-contained instead: they read the same on the dark theme, in the PNG, on GitHub and in a notebook.
- Colours mean the same thing everywhere: the C++ core is blue, Python orange, the C-ABI family violet, WebAssembly and green closure cells green, formats and files aqua, data arrays yellow, regions magenta, refusals and red (fully split) cells red. The values are the `dataviz` reference palette; the ordinal effort ramp on the roadmap map is its sequential blue.
- Corner nodes are solid ink, mid-edge nodes blue, face centres orange, the body centre violet; hidden edges are dashed.
- Text wears ink, never a series colour, except for a one-word annotation that names a coloured element.

## Adding a figure

Write a function returning `Canvas(...).render()` in the right `figures/*.py`, add it to that module's `FIGURES` dict, embed it in a page as `![one-sentence alt text](/diagrams/<name>.svg)`, regenerate, and commit the SVG and PNG. The check fails on an unreferenced or stale figure, so a figure cannot be added and forgotten.
