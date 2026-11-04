# meshio++ for Octave and MATLAB

Octave / MATLAB bindings for [meshio++](https://github.com/loumalouomega/meshioplusplus): read and write the mesh formats used in scientific computing and FEM, and run meshio++'s operations layer (surface extraction, quality, refine, decimate, partition, smooth, slice, isosurface, …) on the results. Every function of the R binding is here as `mio.<name>`.

MIT-licensed, like the rest of meshio++. Developed and tested with GNU Octave (8 or later); the sources use only the MEX API and the MATLAB-compatible language subset, so they also build and run in MATLAB.

## Prerequisite: build and install the C library

```sh
./build/configure.sh --c-api --build
cmake --install build/cpp-release --prefix /your/prefix
```

Make it discoverable with `export PKG_CONFIG_PATH=/your/prefix/lib/pkgconfig` or `export MESHIOPLUSPLUS_HOME=/your/prefix`. At run time the shared library must be on the loader path (`LD_LIBRARY_PATH`, `DYLD_LIBRARY_PATH` or `PATH`).

## Install

Octave:

```sh
tar -C bindings/octave -czf meshioplusplus.tar.gz meshioplusplus
octave --eval "pkg install meshioplusplus.tar.gz"
```

Then `pkg load meshioplusplus` in each session.

MATLAB, or Octave without `pkg`:

```matlab
addpath('bindings/octave/meshioplusplus/inst');
mio_build();          % compiles meshioplusplus_mex once
```

## Use

```matlab
m = mio.read('bracket.msh');
surf = mio.extract_surface(m);
mio.write(surf, 'bracket_surface.vtu');

p = mio.points(surf);            % num_points x 3
t = mio.connectivity(surf, 1);   % num_cells x 3, 1-based
trisurf(t, p(:, 1), p(:, 2), p(:, 3));
```

Arrays have one entity per row, indices are 1-based, data arrays keep their stored class, and meshes are handle objects that free the native mesh when the last reference goes. See `doc/octave.md` (https://loumalouomega.github.io/meshioplusplus/octave) for the full conventions.

## Tests

```sh
cd bindings/octave/meshioplusplus/tests
octave --no-gui --eval "exit(run_tests() > 0)"
```
