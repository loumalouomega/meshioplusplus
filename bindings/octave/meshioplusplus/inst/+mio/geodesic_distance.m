function r = geodesic_distance(mesh, varargin)
%GEODESIC_DISTANCE  Geodesic distance on a triangle surface (heat method).
%   Options: 'sources' (1-based point ids; zeros(0, 1) for none),
%   'source_region' ([]: a region whose points are sources too),
%   'time_factor' (1), 'tolerance' (1e-10) and 'max_iterations' (0 picks
%   one). Returns a struct: mesh (with point data 'geodesic:distance'),
%   num_sources, num_unreached, max_distance, time_step and converged.

    o = mio.internal.opts(varargin, {'sources', zeros(0, 1), 'source_region', [], ...
        'time_factor', 1, 'tolerance', 1e-10, 'max_iterations', 0}, 'mio.geodesic_distance');
    r = meshioplusplus_mex('mio_geodesic_distance', mesh, double(o.sources), ...
        mio.internal.str(o.source_region, 'source_region'), double(o.time_factor), ...
        double(o.tolerance), double(o.max_iterations));
end
