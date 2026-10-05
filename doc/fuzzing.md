# Fuzzing and sanitizers

meshio++ parses files from untrusted sources by design: its readers are reachable from the C ABI, the browser build, a VS Code extension and an MCP server. Two tools keep them honest: an **ASan + UBSan** build of the C++ test suite, run on every pull request, and a **libFuzzer** harness over every native reader, run on a schedule, whose findings become regression inputs that every pull request replays.

## The contract under test

A native reader either returns a mesh or throws `ReadError`. Every entry point enforces the second half: the registry's readers and the Python `_core.*_read` bindings translate any `std::logic_error` or `std::runtime_error` a parser throws (`std::stoll` on a bad token, `.at()` past an end) into a `ReadError` naming the format (`detail/read_guard.hpp`). Anything else is a defect:

- a crash, or a memory error or undefined behaviour reported by ASan/UBSan;
- any other C++ exception leaving a reader;
- a `std::bad_alloc` or an out-of-memory abort, which almost always means a count from a header was trusted before the bytes that should back it were seen;
- a hang.

`detail/parse_guard.hpp` holds the checks a reader makes before trusting its input (`need_tokens`, `checked_count`, `checked_integer`, `file_bytes`), and `NDArray` refuses a shape whose byte size overflows.

## Sanitizer builds

```sh
CC=clang CXX=clang++ build/configure.sh --sanitize address,undefined --tests --c-api --build
ctest --test-dir build/cpp-release-san --output-on-failure
```

`--sanitize` maps to the CMake option `MESHIOPLUSPLUS_SANITIZE` (a `;`-list, for example `address;undefined`), which instruments the core, the C API objects, the GoogleTest binary and the fuzz replay driver. It works with GCC or Clang, and refuses the Python extension: a sanitized `_core` needs the sanitizer runtime preloaded into the interpreter. UBSan is fatal (`-fno-sanitize-recover=undefined`), so a report fails its test. The `cpp-sanitize` CI job runs the whole suite this way with leak checking on; `tools/sanitizers/lsan.supp` is where a third-party leak would be suppressed, and it names none today. Run `ctest` serially: several tests write fixed temporary names, and a parallel run collides on them.

## The fuzz target

```sh
CC=clang CXX=clang++ build/configure.sh --fuzzers --build
python tools/fuzz/seed_corpus.py seeds --replay build/cpp-release-san/meshioplusplus_fuzz_replay
tools/fuzz/run_campaign.sh build/cpp-release-san seeds findings 60
```

`--fuzzers` (CMake `MESHIOPLUSPLUS_BUILD_FUZZERS`, Clang only) builds `meshioplusplus_fuzz_read` with libFuzzer, ASan and UBSan, and compiles the whole core with coverage instrumentation, which is what lets the fuzzer see into the parsers. One binary serves every reader; the format comes from, in order, `MIO_FUZZ_FORMAT`, a `-format=<name>` argument, or the binary's name (`meshioplusplus_fuzz_read_<format>`, so per-format copies work as OSS-Fuzz expects). Each input is written to a per-process scratch file named the way the format is found (its default extension, or a fixed name such as `d3plot` or `z88i1.txt`) and read through `registry_readers()`. Companion-file readers take a deterministic `MIOB` bundle instead (`tests/fuzz/fuzz_bundle.hpp`): XDMF XML/heavy-data sets and ADIOS2 `.bp` directories reconstruct from one bundle input, so findings stay replayable bundle bytes. Library-backed readers (HDF5/netCDF, XDMF bundles, ADIOS2, TecIO, CGNS MLL) run isolated in a forked child so a library abort attributes clearly instead of corrupting the fuzzer's global state; set `MIO_FUZZ_NO_ISOLATE=1` to debug without the fork. Set `MIO_FUZZ_CRASH_ONLY=1` to tolerate every C++ exception while triaging memory errors alone.

`tools/fuzz/seed_corpus.py` builds one seed directory per reader from the fixtures under `tests/python/meshes/` and from a small mesh meshio++ writes in each format, packing XDMF companions (XML, Binary, HDF plus region bundles) and each DOLFINx `.bp` fixture directory into deterministic `MIOB` bundles; `tools/fuzz/dicts/` holds libFuzzer dictionaries for the keyword formats, including XDMF DataItem/reference tokens. `tools/fuzz/run_campaign.sh` runs every format for a fixed time and writes a one-line verdict per format to `summary.txt`; `FORMATS` narrows it, `SHARD`/`NSHARDS` split it.

The native fuzz harness and shell scripts require a POSIX environment. The campaign script supports macOS's Bash 3.2, including formats without dictionaries; Python's shell-script tests skip Windows rather than invoke its WSL launcher. Corpus-packaging tests remain portable.

The ordinary campaign skips `tools/fuzz/not_fuzzed.txt`: Elmer directories and the HDF5/netCDF readers, which need the instrumented dependencies below. XDMF bundles run here (XML, Binary and system-HDF5 HDF paths); the instrumented-library campaign replays the HDF bundles against the sanitized libraries too. ADIOS2, TecIO and CGNS MLL run through the same production-path harness whenever the build enables them (`MESHIOPLUSPLUS_WITH_ADIOS2`/`TECIO`/`CGNSLIB`), with positive bundle seeds and replayable malformed inputs; only Elmer directories remain excluded.

The scheduled `fuzz` workflow runs the campaign weekly in four shards, and uploads any finding with its log.

## Regression inputs

A fixed finding is committed as a minimized input under `tests/fuzz/regressions/<format>/`. With tests enabled, CMake adds one `fuzz_regression_<format>` test per directory, which replays every file through `meshioplusplus_fuzz_replay` -- the same harness with a plain `main()`, so any compiler runs it and a finding can be debugged without libFuzzer. The `cpp-tests` and `cpp-sanitize` jobs therefore replay every past crash on every pull request.

To add one: reproduce with `meshioplusplus_fuzz_read -format=<fmt> <input>`, minimize it with `-minimize_crash=1 -exact_artifact_path=<out>` against the unfixed build, fix the reader, and commit the minimized file under a name that says what it exercises.

## Instrumented library readers

The `libraries` job in the fuzz workflow builds pinned zlib 1.3.1, HDF5 1.14.6 and netCDF 4.9.3 with coverage, ASan and UBSan through `tools/fuzz/build_library_deps.sh`. Source archives are SHA-256 checked, optional networking/plugins and uninstrumented compression dependencies are disabled, and CI checks the actual linked library paths. The prefix records compiler/flags plus the backported patch hash in `instrumentation.txt` and refuses reuse after they change. This is an opt-in Linux development prefix, not a library to load into Python or install system-wide.

HDF5's C library alone is compiled with `-fno-sanitize=function`: its internal `H5I__dec_ref` cleanup callbacks trigger Clang's incompatible-function-pointer check on valid input (observed with Clang 23). ASan and other UBSan checks remain enabled there, and the core, custom mutator, zlib and netCDF retain their full sanitizer flags. This exception is not a general sanitizer-report suppression; revisit it when upgrading HDF5.

The dependency build carries one backported upstream maintenance fix, `tools/fuzz/patches/netcdf-4.9.3-utf8proc-null-guard.patch` (upstream commit `193cf061`, November 2025): the bundled Unicode helper formed a buffer-plus-offset on a NULL base while measuring output length. The check build reported it during valid reads, so the measuring pass now carries NULL through instead. The patch file records its provenance; remove it once the pinned netCDF release includes the fix.

`tools/fuzz/library_formats.txt` names CGNS's HDF5 route, H5M, HMF, MED, MSC Nastran HDF5, VTKHDF, Exodus and XDMF bundles. `meshioplusplus_fuzz_hdf5` reuses the production reader harness and adds bounded mutations through HDF5's own API: numeric values, missing datasets, rank/shape/type changes, attribute replacement and dangling/cyclic soft links. It visits hard-linked objects only and never creates external links. Inputs without the HDF5 magic fall back to ordinary byte mutation, as does classic netCDF; a quarter of HDF5 inputs also use that path. The structure-aware path runs in a forked child so a library abort inside `H5Ovisit`/`H5Dopen` falls back to byte mutation instead of killing the fuzzer. Every library-backed read also runs isolated in a forked child (`MIO_FUZZ_NO_ISOLATE=1` disables it), so an abort inside `H5Fopen`/`nc_open`/ADIOS2/TecIO/cgnslib attributes to the input instead of corrupting the fuzzer's global state, and production readers preflight containers first (`src/cpp/src/formats/library_preflight.hpp`: HDF5/netCDF magic plus the HDF5 superblock version) to refuse obvious corruption as `ReadError` before the library parses it. Findings are ordinary file bytes (or `MIOB` bundle bytes for XDMF/`.bp`) and replay with `meshioplusplus_fuzz_replay`, without a custom decoder or Python fallback.

For a local run, build dependencies into a fresh absolute prefix, point CMake's `CMAKE_PREFIX_PATH`, `HDF5_ROOT` and `ZLIB_ROOT` at it, enable `MESHIOPLUSPLUS_BUILD_FUZZERS`, disable Python and build `meshioplusplus_fuzz_hdf5` and `meshioplusplus_fuzz_replay`. Generate seeds with `seed_corpus.py --library-readers --formats cgns,h5m,hmf,med,nastran_h5,vtkhdf,exodus,xdmf --replay <replay-binary>`. Run the existing campaign with `MIO_FUZZ_LIBRARY_READERS=1`, `FUZZ_BINARY=meshioplusplus_fuzz_hdf5`, `MAX_LEN=262144` and `FORMATS` restricted to those names. For the optional production routes, enable `MESHIOPLUSPLUS_WITH_ADIOS2=ON` (point `CMAKE_PREFIX_PATH` at an ADIOS2 install) and `MESHIOPLUSPLUS_WITH_CGNSLIB=ON` (point `CGNS_ROOT` at a cgnslib install) in the same instrumented build: `read_cgns` then exercises the MLL path first with its raw-HDF5 fallback intact, and `vtx` bundles reconstruct `.bp` directories through ADIOS2's own engine. For TecIO SZL, point `TECIO_ROOT` at a TecIO source tree (Tecplot's distribution or SU2's vendored copy) in the same build: `read_szplt` then decodes through TecIO's own reader into the shared Tecplot zone model, with the existing `.szplt`/`.dat` twins as positive seeds. Each optional route needs only its production decoder plus bundle seeds; no mock decoder counts as coverage. The library job replays seeds before mutation and uploads logs/findings on every run. Use fresh output directories for independent campaigns; an existing finding remains a failure.

The campaign has fixed six meshio defects, each with a native test and a minimized replay input under `tests/fuzz/regressions/`: CGNS coordinate rank/length validation, MED family `NOM` size against `NBR`, VTKHDF connectivity/offset agreement (shape and count), HDF5 link-count/name-length caps in `group_links`/`link_names`, XDMF cell-data size validation in `split_raw_cell_data`, and XDMF HDF datasets larger than their declared `Dimensions`, refused before they are read (a corrupt chunked extent made HDF5 1.10 time out). Native `Cgns.MismatchedCoordinateLengthsThrow`, `Vtkhdf.RefusesConnectivityCountDisagreeingWithIds` and the earlier connectivity-shape probe cover the first three; `tests/python/test_fuzz_tools.py` also generates the shape probe for a supplied `MIO_FUZZ_REPLAY` binary.

Known library limitations remain, distinct from meshio defects: HDF5 aborts on some corrupted dense-link tables (`H5G__dense_iterate`/`H5O__link_reset`), huge object-header allocations, pipeline-copy leaks on failed dataset opens and slow heap deserialization; netCDF-4 inherits the dense-link abort through `nc_open`. These reproduce with the pinned libraries and need upstream fixes; the harness's fork isolation plus the magic/superblock preflight keep them attributable and replayable (a library abort no longer corrupts the fuzzer's global state) while the campaign proceeds. XDMF bundles (XML, Binary, HDF, DataItem references including cycles, mixed topology, sets), ADIOS2 `.bp` bundles (multi-step, mesh-reuse, ghost handling via the existing `test_vtx`/`vtx` fixtures), the CGNS MLL route (ADF containers, NGON_n/NFACE_n via `read_cgns_mll`) and TecIO SZL decoding (shared zone model via `read_szplt` against its `.dat` twins) all run through their production decoders with positive bundle seeds and replayable malformed inputs; `tests/cpp/test_fuzz_bundle.cpp` and `tests/cpp/test_library_preflight.cpp` pin the bundle encoding and the preflight refusals.

The mutator child clears libFuzzer's inherited sanitizer death callback and abort/alarm handlers, so a failed mutation falls back to byte mutation without writing an unrelated empty crash artifact. Reader children retain finding reporting. CI runs `meshioplusplus_fuzz_mutator_probe` through `test_native_mutator_isolation` with injected aborts and ASan errors to verify this distinction; set `MIO_FUZZ_MUTATOR_PROBE` to that executable to run the same tests locally.

## OSS-Fuzz

`tools/fuzz/oss-fuzz/` supplies `Dockerfile`, `build.sh`, seed packaging and `project.yaml` with the maintainer's confirmed Google-account contact. The three upstream integration files carry Apache-2.0 license headers as OSS-Fuzz requires; the meshio++ core and project-owned seeds remain MIT. The image temporarily builds the public `tier-2` branch because the external-engine and seed-generator changes have not yet landed on `master`; switch it to the default branch after they merge, before deleting `tier-2`. It checks out only the Eigen and JSON submodules needed by the core, not the interactive viewer. CMake's `MESHIOPLUSPLUS_FUZZING_ENGINE` accepts OSS-Fuzz's `LIB_FUZZING_ENGINE`; the normal target retains transitive libraries and feature definitions rather than hand-linking object-file globs. With an external engine the project does not add its own sanitizer defaults or coverage flags: OSS-Fuzz's compiler environment controls them. The integration currently disables HDF5/netCDF; the library campaign is separate.

`meshioplusplus_fuzz_seeds` writes small project-generated meshes and validates them through native readers before packaging. `package_seeds.py` combines these with committed regression inputs into per-target `_seed_corpus.zip` archives, rejecting LFS pointers and excluding empty/oversized files. Some read-only or specialised formats have only regression seeds; directory/companion formats are not made meaningful by packaging their files separately.

For container validation, copy `Dockerfile`, `build.sh` and `project.yaml` into `projects/meshioplusplus/` in a checkout of [google/oss-fuzz](https://github.com/google/oss-fuzz). Run the following from that checkout with Docker available. Use `--clean` when changing sanitizers: CMake caches compiler flags, and reusing an address build for an undefined check is not validation of both configurations. Pass the matching sanitizer to `check_build` and `run_fuzzer` too, rather than letting them default to address.

```sh
python3 infra/helper.py build_image meshioplusplus
for sanitizer in address undefined; do
    python3 infra/helper.py build_fuzzers --clean --sanitizer "$sanitizer" meshioplusplus
    python3 infra/helper.py check_build --sanitizer "$sanitizer" meshioplusplus
    python3 infra/helper.py run_fuzzer --sanitizer "$sanitizer" meshioplusplus \
        meshioplusplus_fuzz_read_vtu -- -max_total_time=30
done
python3 infra/helper.py build_fuzzers --clean --sanitizer coverage meshioplusplus
python3 -m zipfile -e \
    build/out/meshioplusplus/meshioplusplus_fuzz_read_vtu_seed_corpus.zip \
    build/corpus/meshioplusplus/meshioplusplus_fuzz_read_vtu
python3 infra/helper.py coverage --no-corpus-download --no-serve \
    --fuzz-target=meshioplusplus_fuzz_read_vtu meshioplusplus
```

The smoke run consumes the packaged positive and regression seeds. For coverage, extract the matching `_seed_corpus.zip` into the target's local `build/corpus/meshioplusplus/<target>/` directory first, or pass that directory through `--corpus-dir`. Inspect the report to confirm that the production reader is reached, not just the harness's input checks. A local external-engine build is useful but is not a substitute for these container checks. No upstream acceptance is implied by the integration files; eligibility and onboarding remain OSS-Fuzz's review decision.
