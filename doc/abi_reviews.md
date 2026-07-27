# ABI additive-change reviews

`tools/check-abi-version.sh` fails a build whose installed headers changed while
`MESHIOPLUSPLUS_ABI_VERSION` did not. This file is how that judgement gets recorded when the
change really was additive.

**One row per ABI version that survived a header change.** Adding a row is a claim that every
listed header change is Tier C by [`doc/abi.md`](abi.md)'s criterion — a *new* inline function,
constexpr variable, type, declaration or header, or an *appended* enumerator — and therefore
cannot affect a consumer that was already compiled. It is not a rubber stamp: the whole point
of the ABI number is that a consumer may trust it instead of re-pinning the release version,
so a wrong entry here ships silent memory corruption.

If the change edited the *body* of an existing inline function or template, or moved a data
member, it is Tier A/B and belongs in a version bump, not in this table.

| ABI | releases | headers changed | why it is additive |
| --- | --- | --- | --- |
| 3 | v9.3.0 | `formats/exodus.hpp` | Adds one new `inline constexpr const char* kExodusAttributePrefix`. No existing declaration touched, `ExodusInfo` unchanged. A consumer compiled against v9.2.0 headers has an identical translation unit either way — it simply does not have the name. |
| 3 | v9.4.0 | `abi_version.hpp`, `detail/abi_version_check.hpp`, `mesh.hpp` | All new: the ABI machinery itself. `mesh.hpp` gains only an `#include`; no type's definition changes. See the note in [`doc/abi.md`](abi.md#the-one-false-positive) about the one-way incompatibility this introduction creates (v9.4.0 headers reference a sentinel a ≤v9.3.0 library does not define), which fails closed at link time. |
