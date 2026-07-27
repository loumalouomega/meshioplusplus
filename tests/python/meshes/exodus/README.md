# Exodus reference fixtures

## `DCBmodel_PD_solid.e`

A real peridynamics output file: a 2-D double-cantilever-beam model discretized
as 504 one-node `SPHERE` particles in four element blocks, with two node sets,
nine nodal fields and ten time steps.

**This is not a synthetic fixture, and that is the point.** It is the file from
[VSCode-MDPA-Preview#63](https://github.com/loumalouomega/VSCode-MDPA-Preview/issues/63),
where meshio++ failed with `Exodus: unknown element type SPHERE` on a type it had
mapped correctly all along. The cause was in the *encoding*, not the type table:
[NetCDF.jl](https://github.com/JuliaGeo/NetCDF.jl), which PeriLab writes Exodus
with, counts the C string's terminating NUL as part of the `elem_type`
attribute's length, so the value arrives as the 7 characters `"SPHERE\0"`.
`netCDF4` strips that on the way in, which is exactly why no amount of
Python-side testing — round-trip or otherwise — could have caught it, and why
the failure only ever showed up in WASM, where there is no Python to fall back
to. Nothing meshio++ can generate itself reproduces the property; only a file
written by that toolchain does.

The hand-authored `write_sphere_fixture` in `tests/python/exodus_fixture.py`
reproduces the *shape* (and can vary it, which this file cannot), but it has to
patch the classic-format bytes by hand to get the NUL in. This file is the
ground truth those constructions are checked against.

| | |
|---|---|
| **Source** | [`PeriHub/PeriLab.jl`](https://github.com/PeriHub/PeriLab.jl), `test/fullscale_tests/test_DCB/Reference/DCBmodel_PD_solid.e` |
| **Retrieved** | 2026-07-27, from `main` |
| **SHA-256** | `c0be64f9949a918126b16934837e2ca0e0de2f62fc9ef7d920e5b2c8ae9703d3` |
| **Copyright** | 2023 Christian Willberg, Jan-Timo Hesse (DLR) |
| **Licence** | BSD-3-Clause — see `LICENSE.BSD-3-Clause` and the `.license` sidecar |
| **Modified** | No; committed byte-for-byte as retrieved |

BSD-3-Clause is permissive and imposes no licence change on this MIT repository,
but its first condition does require the copyright notice, the conditions and
the disclaimer to travel with the file — hence `LICENSE.BSD-3-Clause` here and
the upstream [REUSE](https://reuse.software/) `.license` sidecar kept verbatim
beside the data. Do not move the file without them.
