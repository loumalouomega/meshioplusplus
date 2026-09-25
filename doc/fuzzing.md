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

`--fuzzers` (CMake `MESHIOPLUSPLUS_BUILD_FUZZERS`, Clang only) builds `meshioplusplus_fuzz_read` with libFuzzer, ASan and UBSan, and compiles the whole core with coverage instrumentation, which is what lets the fuzzer see into the parsers. One binary serves every reader; the format comes from, in order, `MIO_FUZZ_FORMAT`, a `-format=<name>` argument, or the binary's name (`meshioplusplus_fuzz_read_<format>`, so per-format copies work as OSS-Fuzz expects). Each input is written to a per-process scratch file named the way the format is found (its default extension, or a fixed name such as `d3plot` or `z88i1.txt`) and read through `registry_readers()`. Set `MIO_FUZZ_CRASH_ONLY=1` to tolerate every C++ exception while triaging memory errors alone.

`tools/fuzz/seed_corpus.py` builds one seed directory per reader from the fixtures under `tests/python/meshes/` and from a small mesh meshio++ writes in each format; `tools/fuzz/dicts/` holds libFuzzer dictionaries for the keyword formats. `tools/fuzz/run_campaign.sh` runs every format for a fixed time and writes a one-line verdict per format to `summary.txt`; `FORMATS` narrows it, `SHARD`/`NSHARDS` split it.

Readers listed in `tools/fuzz/not_fuzzed.txt` are skipped: `elmer` reads a directory, and the HDF5, netCDF, ADIOS2 and TecIO formats hand their bytes to a third-party library that is neither instrumented here nor this project's to fix.

The scheduled `fuzz` workflow runs the campaign weekly in four shards, and uploads any finding with its log.

## Regression inputs

A fixed finding is committed as a minimized input under `tests/fuzz/regressions/<format>/`. With tests enabled, CMake adds one `fuzz_regression_<format>` test per directory, which replays every file through `meshioplusplus_fuzz_replay` -- the same harness with a plain `main()`, so any compiler runs it and a finding can be debugged without libFuzzer. The `cpp-tests` and `cpp-sanitize` jobs therefore replay every past crash on every pull request.

To add one: reproduce with `meshioplusplus_fuzz_read -format=<fmt> <input>`, minimize it with `-minimize_crash=1 -exact_artifact_path=<out>` against the unfixed build, fix the reader, and commit the minimized file under a name that says what it exercises.

## OSS-Fuzz

`tools/fuzz/oss-fuzz/` holds a draft `project.yaml` and `build.sh` for an OSS-Fuzz submission; the submission itself is an open item on the [roadmap](./roadmap.md).
