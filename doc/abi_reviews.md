# ABI additive-change reviews

`tools/check-abi-version.sh` fails a build whose installed headers changed while
`MESHIOPLUSPLUS_ABI_VERSION` did not. This file is how that judgement gets recorded when the
change really was additive.

**One row per release that changed an installed header while holding the ABI version.** Adding a
row is a claim that every listed header change is Tier C by [`doc/abi.md`](abi.md)'s criterion —
a *new* inline function, constexpr variable, type, declaration or header, or an *appended*
enumerator — and therefore cannot affect a consumer that was already compiled. It is not a
rubber stamp: the whole point of the ABI number is that a consumer may trust it instead of
re-pinning the release version, so a wrong entry here ships silent memory corruption.

**The `releases` and `headers changed` columns are read by the gate, not just by people.**
`tools/check-abi-version.sh` looks for a row matching the ABI version *and* the release version
in `CMakeLists.txt`, then requires every header that actually changed to be named in it — paths
relative to `src/cpp/include/meshioplusplus/`, as written below. A header you changed but did
not list is a red build naming it.

That specificity is the fix for a real defect: through v9.4.0 the gate matched any row whose
first column was the current ABI number, so once *any* ABI-3 row existed, every later header
change matched it and passed — reporting "records the additive review" about a review of a
different release entirely. A row for one release says nothing about the next one.

If the change edited the *body* of an existing inline function or template, or moved a data
member, it is Tier A/B and belongs in a version bump, not in this table.

| ABI | releases | headers changed | why it is additive |
| --- | --- | --- | --- |
| 3 | v9.3.0 | `formats/exodus.hpp` | Adds one new `inline constexpr const char* kExodusAttributePrefix`. No existing declaration touched, `ExodusInfo` unchanged. A consumer compiled against v9.2.0 headers has an identical translation unit either way — it simply does not have the name. |
| 3 | v9.4.0 | `abi_version.hpp`, `detail/abi_version_check.hpp`, `mesh.hpp` | All new: the ABI machinery itself. `mesh.hpp` gains only an `#include`; no type's definition changes. See the note in [`doc/abi.md`](abi.md#the-one-false-positive) about the one-way incompatibility this introduction creates (v9.4.0 headers reference a sentinel a ≤v9.3.0 library does not define), which fails closed at link time. |
| 4 | v9.7.0 | `formats/gmsh.hpp`, `version.hpp` | `version.hpp` is the routine release bump (macro values only, no declaration). In `formats/gmsh.hpp`, all new: a `GmshInfo` struct (the `MedInfo`/`ExodusInfo` side-channel pattern) and two overloads that take it, `read_gmsh(path, GmshInfo&, ReadOptions)` and `write_gmsh41(path, mesh, binary, const GmshInfo&)`. No existing declaration, default argument or type definition is touched, and neither existing function's mangled name changes — a consumer compiled against v9.6.0 headers links against exactly the same symbols. The one source-level (not ABI) consequence is that `&read_gmsh` and `&write_gmsh41` are now ambiguous as bare function pointers, which is a **compile** error at the call site, never silent; `registry.cpp` wraps both in lambdas for that reason. |
