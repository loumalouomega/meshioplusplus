# Dataset dashboard

The [dataset manager](https://loumalouomega.github.io/meshioplusplus/viewer/dataset.html) — the browser viewer's second page — curates the [dataset manifests](./datasets) a training run is driven off. Since v10.22.0 it has two depths: an **overview** of every manifest in the picked directory, and the per-manifest **curation view** that page always had. Everything on this page runs client-side against the viewer's own WASM build; no file is uploaded anywhere. Launching and monitoring PhysicsNeMo training from here needs a local companion process and is the next item of [the roadmap](./roadmap) — it is not built yet.

![The dataset dashboard: one card per manifest in the workspace, with split balance, tags, health badges and a preview thumbnail](/viewer/dataset-dashboard.png)

## Two depths of one page

Picking a case directory (the File System Access picker in Chromium, a `webkitdirectory` input elsewhere — see [the per-browser table](./datasets#curating-in-the-browser)) lands on the overview. Every `*.json` file at the workspace root is parsed as a manifest and shown as a card; **Open** drills into that manifest's curation view (entries, detail editor, preview, scan), and **← All datasets** returns to the grid. Unsaved edits survive the round trip — the card of the manifest being edited is flagged, its counts and split balance follow the in-memory document, and opening it again does not re-read the file underneath your changes. `Scan all` refuses to run while a manifest has unsaved edits, since it opens every manifest in turn.

The page's state hook, `window.__datasetState`, records which depth is showing (`view`), the cards (`manifests`) and the last diff (`diff`), which is what the Playwright suite asserts on.

## Overview cards

Each card shows the manifest's name (or path), its description, the entry count and the file's last-modified time, a proportional **split-balance bar** with a colour-and-label legend (split colours are assigned in a fixed order by sorted split name across the whole workspace, so `train` is the same colour on every card; unassigned entries are the grey tail), the union of its tags and groups as chips, its health badges (below), and a thumbnail. Cards sort by name, entry count, last-modified time or health (worst first).

A manifest the strict parser rejects — an unknown key, a duplicate id, invalid JSON — is still a card, flagged with the parser's message, rather than a silently missing one: a broken manifest is exactly what an overview should surface. **New manifest** starts an empty `dataset.json` in the curation view; it gets a card of its own until it is saved.

Thumbnails are captured from the preview render — the first time an entry of a manifest is previewed, and during **Scan all** when *thumbnails* is ticked — and live only for the session; a card only ever lacks a thumbnail, never shows a stale one.

## Health summaries

**Scan** (in a manifest's curation view, or the card's own Scan button) stages every entry once, one at a time, and summarizes its data arrays and quality metrics at step 0. The per-entry scan records the step count, NaN/Inf counts over the *data* arrays, the inverted and degenerate cell counts, the worst scaled Jacobian and the list of data arrays present. Scans are cached (in IndexedDB, keyed by the manifest, the entry and the newest modification time of its source files), so re-opening a manifest shows its health without re-staging anything and a changed source file invalidates only its own entry.

One rule is load-bearing, and it has a single home (`src/viewer/src/dataset/health.ts`): a `quality:*` metric's NaN means *"does not apply to this cell type"* by design — `compute_quality`'s own convention — so quality rows are **excluded** from the NaN/Inf counts; a quality row flags red only on actual inverted or degenerate cells or a negative worst scaled Jacobian.

The manifest-level summary aggregates the scans into badges shown on the card and in the curation view's **Health** section:

| badge | meaning |
|---|---|
| `scanned k/n` | how many entries have a scan (a partial scan is a warning) |
| `unassigned k` | entries with no split |
| `NaN/Inf n`, `inverted n`, `degenerate n` | totals over the scanned entries — a bad dataset is visible before a run wastes GPU time on it |
| `min SJ x` | the worst scaled Jacobian (red below 0, amber below 0.2) |
| `fields missing in k` | entries lacking a data array that every other scanned entry carries — the Health section lists which arrays; a mismatched field set is what makes `feature_matrix` fail at training time |
| `unreadable n` | entries that could not be staged or summarized at all |
| `healthy` | a fully scanned manifest with none of the above |

Everything here is computed in the browser from the same summaries the per-entry table shows. A server-side producer of the same summary — for manifests too large to stage in a browser, and for the companion process the training features need — is the roadmap's next step and will fill the same `ManifestHealth` shape with `producer: "server"`.

## Diffing manifests

**Diff…** compares two manifests side by side. Either side can be the open manifest **on disk**, its **current edits**, any manifest at the workspace root, or **pasted JSON**. The result is structural first — document fields that changed, entries added, entries removed, and for each entry present in both the fields that differ (split, tags, group, notes, metadata, source; metadata is compared by value, so key order never reads as a change) — followed by a line diff of the two serializations with unchanged context collapsed. Because the manifest serialization is stable (`indent=2`, fixed key order, trailing newline), the line diff of a hand edit is as clean as `git diff`'s.

Git history is out of a browser page's reach: to compare against a revision, paste the output of `git show REV:dataset_manifest.json` as one side.

## Limitations

- Only `*.json` files at the workspace root are treated as manifests: the manifest is the resolution anchor for the relative sources inside it, and the page saves back at the root.
- Thumbnails are session-only and cost one preview render per manifest; untick *thumbnails* before **Scan all** on a large workspace.
- A scan stages exactly one entry at a time (the WASM filesystem never reclaims memory), so **Scan all** over many large manifests is serial and slow the first time; the scan cache makes the next visit instant.
- Training launch, monitoring, run history, prediction preview and notifications are not built yet — see [the roadmap](./roadmap).
