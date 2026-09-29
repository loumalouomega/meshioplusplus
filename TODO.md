# TODO — roadmap §4 PR1 (v16.27.0, ABI 20), branch `feat/core-parity-1`

Plan: `~/.claude/plans/start-with-tier-4-peppy-naur.md`. Memory: `core-parity-task-state.md`. Delete this file before merging.

## State
Code is done and committed: MDPA `MdpaInfo` side channel + `mio_format_info` (C/Fortran/Julia/R/WASM/pybind), Side regions surviving operations (`remap_region`, `_side_carry.py`), XDMF `<Set>`s and the VTU/VTP `<FieldData>` region convention. Master (v16.26.0 fuzz fixes) is merged; this is v16.27.0. Version files, CHANGELOG, `doc/abi.md` row 20, roadmap, `BASELINE_HASHES` and the single header are refreshed. ABI gate OK.

## To do
1. `git push -u origin feat/core-parity-1`, then open the PR against master. No "Generated with Claude Code" / `Co-Authored-By` lines (user rule). Mention the unverified items below in the body.
2. Re-run on the v16.27.0 build (a final run was interrupted): gtest, Fortran (`~/.cache/mio-scratch-p4/fortran.sh`), full pytest (`~/.cache/mio-scratch-p4/verify.sh` / `final.sh`). gtest passed; NATIVE/KRATOS passed except two HDF5 `XdmfTimeSeries` baseline failures.
3. Determinism sweep vs a master baseline: `tools/bench_ops.sh` with `--hash` (see AGENTS.md). Expect no digest change (no golden fixture carries regions).
4. Docs build in `doc/` (`npm install`, `npm run docs:build`) and check links by hand (`ignoreDeadLinks` is on).
5. Lint: `pre-commit run -a` (or `isort -p meshioplusplus`, black, flake8, clang-format on touched hunks). `tools/include-cleanup.sh --check`.
6. Not verified locally: WASM (`bindings/wasm/js_bindings.cpp`, `tests/wasm/smoke.mjs`; no Emscripten), Julia (parse-checked only), R glue (syntax-checked only).
7. Example notebook demonstrating Side-region carry and the new formats (`example/python/`, PyVista off-screen + matplotlib fallback, executed and committed) — not written.
8. Fuzz seeds for the new parsers (MDPA info path, XDMF `<Set>`, region `<FieldData>`) and `fuzz_read.cpp` calling the `MdpaInfo` overload — not done.
9. MCP: no public Python API signature changed, so `_tools.py`/`_server.py` only got docstring edits; confirm `tests/python/test_mcp.py` still passes.

## Later PRs (stacked, same ABI 20; do not tag until the stack lands)
- v16.28.0: Gmsh `$Periodic` + 4.0 in C++, VTK family (appended/multi-piece VTP, VTS/VTR/VTI appended, legacy structured), Netgen extras, XDMF 2 + `Reference="XML"`.
- v16.29.0: MED multi-mesh/profiles, Exodus sets + series writer.
- v16.30.0: pipeline report JSON, PCD `binary_compressed`, glTF options, sets↔data (all additive).
