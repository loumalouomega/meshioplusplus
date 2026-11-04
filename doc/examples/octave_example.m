% Minimal meshio++ Octave / MATLAB consumer: build a tet mesh, write it, read
% it back.
%
% Run with the binding on the path and the C library on the loader path, e.g.
%   LD_LIBRARY_PATH=/opt/meshioplusplus/lib octave --eval ...
%     "addpath('bindings/octave/meshioplusplus/inst'); run('doc/examples/octave_example.m')"

fprintf('meshio++ %s (backend: %s)\n', mio.version(), mio.mesh_backend());

% One point per row (num_points x dim) and one cell per row
% (num_cells x nodes_per_cell), the MATLAB convention: trisurf/patch take
% these directly. Connectivity is 1-based here.
points = [0.0, 0.0, 0.0; 1.1, 0.2, 0.3; 0.4, 1.2, 0.5; 0.6, 0.7, 1.3; 1.4, 1.5, 1.6];
conn = [1, 2, 3, 4; 2, 3, 4, 5];

path = [tempname() '.vtu'];

m = mio.mesh();
mio.set_points(m, points);
mio.add_cell_block(m, 'tetra', conn);
mio.add_point_data(m, 'id', int32((1:5)'));
mio.write(m, path);
mio.release(m);

r = mio.read(path);
assert(mio.num_points(r) == 5);
assert(isequal(mio.points(r), points));
assert(strcmp(mio.cell_block_type(r, 1), 'tetra'));
assert(isequal(mio.connectivity(r, 1), conn));
% The _raw reader is 0-based -- the ABI's own numbering. It still copies.
assert(isequal(mio.connectivity_raw(r, 1), conn - 1));
% Data arrays keep their stored class.
assert(isa(mio.point_data(r, 'id'), 'int32'));
clear r; % the handle frees the native mesh when its last reference goes
delete(path);

fprintf('octave_example.m: OK\n');
