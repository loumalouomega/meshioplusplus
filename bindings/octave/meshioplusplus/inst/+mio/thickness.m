function r = thickness(mesh, varargin)
%THICKNESS  Wall thickness or clearance per facet.
%   R = mio.thickness(MESH, ...) casts a ray from each 2-D cell's centroid
%   along its normal (roadmap 3.1.4): 'inward' stops at the far side of the
%   wall, 'outward' at the next wall facing the ray, and 'auto' is inward for
%   one surface, outward against 'other'. The surface must be consistently
%   oriented with outward normals. Options:
%     'other'        the surface to reach, or [] for MESH itself ([])
%     'direction'    'auto', 'inward' or 'outward' ('auto')
%     'max_distance' the longest ray, 0 for the bounding-box diagonal (0)
%     'threshold'    count the hits below it, 0 for none (0)
%     'region'       named cell region measured ('')
%     'region_other' named cell region reached ('')
%   R is a struct: mesh (with the 'thickness:distance' and
%   'thickness:hit_cell' cell data), direction, num_cells, num_hit,
%   num_missed, num_below, min, max, mean, min_cell (1-based, 0 when nothing
%   was hit), min_point (1 x 3), max_distance and warnings (n x 1 cellstr).

    o = mio.internal.opts(varargin, {'other', [], 'direction', 'auto', 'max_distance', 0, ...
        'threshold', 0, 'region', '', 'region_other', ''}, 'mio.thickness');
    r = meshioplusplus_mex('mio_thickness', mesh, o.other, ...
        mio.internal.str(o.direction, 'direction'), double(o.max_distance), ...
        double(o.threshold), mio.internal.str(o.region, 'region'), ...
        mio.internal.str(o.region_other, 'region_other'));
    r.warnings = split_warnings(r.warnings);
end

function w = split_warnings(text)
    if isempty(text)
        w = cell(0, 1);
    else
        w = strsplit(text, char(10));
        w = w(:);
    end
end
