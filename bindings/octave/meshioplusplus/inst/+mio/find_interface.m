function r = find_interface(mesh, region_a, region_b, varargin)
%FIND_INTERFACE  Matching boundary facets of two Cell regions.
%   R = mio.find_interface(MESH, REGION_A, REGION_B, ...) with options
%   'mesh_b' ([]: both regions are in MESH), 'mode' ('conforming' or
%   'proximity'), 'master' ('a' or 'b': the returned facet mesh and
%   signed-gap normal), 'gap_tolerance' (0), 'angle_tolerance' (30) and
%   'overlap_tolerance' (0). Returns a struct: mesh, num_pairs, area,
%   max_gap, unmatched_a, unmatched_b, side_a and side_b -- the side
%   matrices are (num_pairs x 2), one (cell, facet) pair per row: 1-based
%   cell ids, 0-based local-facet ordinals.

    o = mio.internal.opts(varargin, {'mesh_b', [], 'mode', 'conforming', 'master', 'a', ...
        'gap_tolerance', 0, 'angle_tolerance', 30, 'overlap_tolerance', 0}, ...
        'mio.find_interface');
    r = meshioplusplus_mex('mio_find_interface', mesh, mio.internal.str(region_a, 'region_a'), ...
        mio.internal.str(region_b, 'region_b'), o.mesh_b, mio.internal.str(o.mode, 'mode'), ...
        mio.internal.str(o.master, 'master'), double(o.gap_tolerance), ...
        double(o.angle_tolerance), double(o.overlap_tolerance));
end
