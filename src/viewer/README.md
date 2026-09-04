# meshio++ browser viewer

A client-side mesh viewer and format converter, built on the
[`@meshioplusplus/wasm`](../wasm) package and [vtk.js](https://kitware.github.io/vtk-js/).

Live at **<https://loumalouomega.github.io/meshioplusplus/viewer/>**.

Everything runs in the browser: no server, no upload. A file you open never
leaves your machine.

## How it works

```
file  ──►  Web Worker ──────────────────────────────────►  main thread
           meshio++ WASM                                   vtk.js
           readMetadata()  → the info panel                XMLPolyDataReader
           convertSurface() → VTP  ── transferred ──►      mapper + scalar bar
           convert()        → any writable format
```

Two things drive that shape.

**vtk.js has no unstructured-grid model.** `vtkPolyData` is its only mesh type,
and `XMLPolyDataReader` its only mesh reader — there is no
`XMLUnstructuredGridReader`. So VTP, not VTU, is the interchange format here,
and a volume mesh is displayed by its boundary surface. That is not a
limitation the viewer chose; it is what a surface renderer can draw.

**The JS mesh representation is flat and lossy.** `readMesh` hands back flat
typed arrays that cannot carry a multi-component (vector or tensor) array. So
the viewer never round-trips a mesh through JavaScript: `convertSurface` reads,
extracts the boundary, linearizes and writes VTP entirely inside C++, and
`convert` does the same for the download path. A `readMesh` → `extractSkin` →
`writeMesh` pipeline would have silently dropped every vector field on the way
to the renderer.

## Development

```sh
npm install
npm run dev          # http://localhost:5173/meshioplusplus/viewer/
```

The viewer depends on `file:../wasm`, the package built from this repository.
**npm *copies* a `file:` dependency rather than symlinking it — and a later
`npm install` does not refresh that copy**, because nothing in `package.json`
changed so npm skips it. After touching anything under `bindings/wasm/` or
`src/cpp/`, rebuild the package and force the copy:

```sh
source /path/to/emsdk/emsdk_env.sh
./build/configure-wasm.sh --build                      # repository root
cd src/viewer && rm -rf node_modules/@meshioplusplus && npm install
```

`vite.config.mjs` checks this at `buildStart` and fails with those exact
instructions. It compares two things, because there are two failure modes: the
wrapper's exports (a *missing* binding shows up at runtime as
`m.convertSurface is not a function`) and the `.wasm` binary itself (a
*changed* binding shows up as nothing at all — the app works, quietly using
the old behaviour). The second one is the expensive one.

CI never hits this: it builds the package with emsdk and then runs `npm ci`,
which wipes `node_modules` first.

## Pages

The web build ships two pages: `index.html` (the viewer) and `dataset.html`
(the dataset manager, `src/dataset/`). The latter has two depths — an
overview of every manifest in the picked directory (cards, health, diff) and
one manifest's curation view — documented in `doc/dashboard.md`; its
dashboard-only CSS lives in `src/dataset/style.css`, never in the shared
`src/style.css`, so the embed build's bytes cannot drift.

## Builds

| script | output | used by |
| --- | --- | --- |
| `npm run build:web` | `dist-web/` | the public demo, deployed to `/meshioplusplus/viewer/` by `docs.yml` |
| `npm run build:embed` | `dist-embed/index.html` | a single self-contained file shipped in the Python wheel |

The embedded build backs `meshioplusplus.view(backend="browser")`. There Python
already has the C++ core and has written the VTP itself, so the page only has to
render: the WASM client is aliased to a stub, and the worker and the ~2 MB
`.wasm` never enter the bundle. Being a single file is also what lets it open
over `file://` with no local HTTP server.

Regenerate the committed wheel asset with:

```sh
npm run build:embed
cp dist-embed/index.html ../python/meshioplusplus/_viewer_assets/viewer.html
```

## Tests

```sh
npx playwright install chromium
npm run test:e2e
```

The smoke tests drive the real app against the real WASM package and assert on
`window.__viewerState` rather than on pixels, so they need no GPU. That is
deliberate: nothing else in CI would notice if the published package's API
drifted away from what the viewer calls.
