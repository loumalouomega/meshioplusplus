function r = deviation_regions(mesh, tolerance, varargin)
%DEVIATION_REGIONS  Excess and missing-material regions beyond a tolerance.
%   R = mio.deviation_regions(MESH, TOLERANCE, ...) finds the connected
%   regions (cells joined across shared edges) of a measured surface whose
%   deviation, the mean of each cell's corners' 'field' values, lies beyond
%   +TOLERANCE (excess) or -TOLERANCE (missing material) (roadmap 7.4.3), and
%   writes them as named cell regions with the 'deviation:region' cell data.
%   Options:
%     'field'    the signed deviation point data ('deviation:distance')
%     'prefix'   region names are <prefix>excess_<k>, <prefix>missing_<k>
%                ('deviation:')
%     'min_area' drop regions smaller than this area (0)
%   R is a struct: mesh, regions (n x 1 struct array: name, sign, num_cells,
%   area, extreme, mean, volume, centroid 1 x 3, bounds 1 x 6; largest
%   |volume| first), num_excess, num_missing, num_dropped and warnings.

    o = mio.internal.opts(varargin, {'field', 'deviation:distance', ...
        'prefix', 'deviation:', 'min_area', 0}, 'mio.deviation_regions');
    r = meshioplusplus_mex('mio_deviation_regions', mesh, double(tolerance), ...
        mio.internal.str(o.field, 'field'), mio.internal.str(o.prefix, 'prefix'), ...
        double(o.min_area));
    if isempty(r.warnings)
        r.warnings = cell(0, 1);
    else
        w = strsplit(r.warnings, char(10));
        r.warnings = w(:);
    end
end
