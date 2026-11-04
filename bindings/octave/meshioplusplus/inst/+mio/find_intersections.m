function r = find_intersections(mesh, varargin)
%FIND_INTERSECTIONS  Where a surface intersects itself, or two surfaces meet.
%   R = mio.find_intersections(MESH, ...) reduces the surfaces to triangles
%   and decides every candidate pair with an exact triangle-triangle test
%   (roadmap 3.1.2). On one surface, triangles sharing corners count only
%   when they meet beyond them; against 'other', any contact counts, touching
%   included. Run mio.extract_surface on a volume mesh first. Options:
%     'other'        a second surface, or [] to check MESH against itself ([])
%     'region'       a named cell region of MESH to restrict to ('')
%     'region_other' the same for 'other' ('')
%     'segments'     build the line mesh of the intersection segments (true)
%     'max_pairs'    keep at most this many pairs, 0 for all (0)
%   R is a struct: mesh (with the 'intersections:count' cell data), other
%   (likewise, or []), segments (a line mesh, or []), pairs (num_pairs x 2,
%   1-based global cell ids, one pair per row), num_pairs,
%   num_triangle_pairs, num_coplanar, num_segments, num_cells_a, num_cells_b,
%   truncated, clean and warnings (n x 1 cellstr).

    o = mio.internal.opts(varargin, {'other', [], 'region', '', 'region_other', '', ...
        'segments', true, 'max_pairs', 0}, 'mio.find_intersections');
    r = meshioplusplus_mex('mio_find_intersections', mesh, o.other, ...
        mio.internal.str(o.region, 'region'), mio.internal.str(o.region_other, 'region_other'), ...
        logical(o.segments), double(o.max_pairs));
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
