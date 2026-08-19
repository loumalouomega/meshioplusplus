# Icon sources

TikZ-drawn toolbar icons for the browser viewer, sharing their visual language with the sibling [CAD-Preview](https://github.com/loumalouomega) and MDPA-Preview extensions — thirteen of the eighteen are adapted from MDPA-Preview's set so the three projects look like one family.

All sources share `preamble.tex`: `line width=1.3pt`, round caps and joins, a canvas of roughly -13..13 mm.

## Pipeline

```
tikz-ui/<id>.tex → pdflatex → dvisvgm → svg-ui/<id>.svg
                 → build-icons.mjs → ../../src/viewer/src/ui/icons.ts
```

The same `pdflatex` + `dvisvgm` pair [`logo/build.sh`](../logo/build.sh) already uses. Note that the sibling projects use `pdftocairo` instead, so the post-processing regexes are **not** interchangeable — dvisvgm emits single-quoted attributes, hex colours, `pt` dimensions, and wraps everything in `<g id='page1'>`.

```bash
cd doc/icons
make ui     # tex → svg   (needs pdflatex + dvisvgm)
make ts     # svg → ../../src/viewer/src/ui/icons.ts   (pure Node)
make clean
node check-icons.mjs   # what CI runs; needs neither TeX nor a build
```

`svg-ui/*.svg` is committed so the TypeScript can be regenerated without a TeX install — only changing a drawing needs one.

## Sentinel colours

Icons are drawn in three deliberately odd hex colours that `build-icons.mjs` maps to `currentColor`:

| source colour | becomes |
| --- | --- |
| `iconfg` `#102030` | `currentColor` |
| `iconfg60` `#405060` | `currentColor` at 60% opacity |
| `iconfg30` `#708090` | `currentColor` at 30% opacity |

Drawing in plain black would have meant guessing how dvisvgm encodes it — `#000`? `#000000`? omitted, relying on the SVG default? — so sentinels turn the mapping into an exact string match. They also give shaded faces a way to stay *relatively* shaded while still following the surrounding text colour.

## The `id` strip

`build-icons.mjs` removes every `id` attribute. dvisvgm wraps its output in `<g id='page1'>`, and eighteen icons inlined into one document would mean eighteen elements with `id="page1"` — invalid DOM, and a `getElementById` collision in an app whose entire Playwright hook is `getElementById`.

## Adding an icon

Create `tikz-ui/<newId>.tex` (see `_template.tex.example`), run `make ts`, then reference `ICONS.newId`. Because `IconId` is a union, `tsc --noEmit` checks every reference — a typo is a compile error, not a blank square.

**Never hand-edit `src/viewer/src/ui/icons.ts`**; `make ts` rewrites it wholesale.
