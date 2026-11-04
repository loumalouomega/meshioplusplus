function m = hex_block(n, spacing)
%HEX_BLOCK  An n x n x n hexahedron block on the lattice [0, n*spacing]^3.
%   M = hex_block(N) or hex_block(N, SPACING) builds the small, easy-to-read
%   fixture several operation demos use (convert_cells, refine, smooth,
%   interpolate, slice, isosurface), mirroring the C++, R and Julia
%   notebooks' hex_block. Notebook helper, not part of the package API.

    if nargin < 2
        spacing = 1.0;
    end
    % Point (i, j, k), each 0..n, has the 1-based index i*(n+1)^2 + j*(n+1) + k + 1.
    pid = @(i, j, k) i * (n + 1)^2 + j * (n + 1) + k + 1;
    [k, j, i] = ndgrid(0:n, 0:n, 0:n);
    pts = [i(:), j(:), k(:)] * spacing;         % num_points x 3

    [k, j, i] = ndgrid(0:n - 1, 0:n - 1, 0:n - 1);
    i = i(:);
    j = j(:);
    k = k(:);
    conn = [pid(i, j, k), pid(i + 1, j, k), pid(i + 1, j + 1, k), pid(i, j + 1, k), ...
            pid(i, j, k + 1), pid(i + 1, j, k + 1), pid(i + 1, j + 1, k + 1), pid(i, j + 1, k + 1)];

    m = mio.mesh();
    mio.set_points(m, pts);
    mio.add_cell_block(m, 'hexahedron', conn);
end
