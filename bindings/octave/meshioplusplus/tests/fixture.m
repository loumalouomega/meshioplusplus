function m = fixture()
%FIXTURE  The shared test mesh: 5 points x 3 dims, 2 tetra x 4 nodes,
%   3-component vector data, coordinates asymmetric in every axis -- the
%   same fixture as the Fortran, Julia and R suites. Nothing is square, so a
%   transposed mapping or a missed 1-based shift cannot cancel out.

    m = mio.mesh();
    mio.set_points(m, fixture_points());
    mio.add_cell_block(m, 'tetra', fixture_conn());
    mio.add_point_data(m, 'displacement', fixture_vec());
    mio.add_point_data(m, 'temperature', [1; 2; 3; 4; 5]);
    mio.append_cell_data(m, 'material', [7; 9]);
end
