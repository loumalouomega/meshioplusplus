function m = cube_surface()
%CUBE_SURFACE  The unit cube as a closed, outward-wound triangle surface.
    m = mio.mesh();
    mio.set_points(m, [0 0 0; 1 0 0; 1 1 0; 0 1 0; 0 0 1; 1 0 1; 1 1 1; 0 1 1]);
    mio.add_cell_block(m, 'triangle', [1 3 2; 1 4 3; 5 6 7; 5 7 8; 1 2 6; 1 6 5; ...
                                       2 3 7; 2 7 6; 3 4 8; 3 8 7; 4 1 5; 4 5 8]);
end
