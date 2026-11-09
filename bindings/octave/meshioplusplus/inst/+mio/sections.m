function m = sections(mesh, direction, offsets, varargin)
%SECTIONS  Chained sections of a surface at parallel planes.
%   M = mio.sections(MESH, DIRECTION, OFFSETS) slices MESH at the planes
%   n . x = OFFSETS(k) (n the unit DIRECTION), walks each section into
%   polylines and merges them into one line mesh with the section:offset,
%   section:chain and section:closed cell data (roadmap 7.6.1). Options:
%   'min_length' (0) and 'record_parent_ids' (false).

    o = mio.internal.opts(varargin, {'min_length', 0, 'record_parent_ids', false}, ...
        'mio.sections');
    m = meshioplusplus_mex('mio_sections', mesh, double(direction), double(offsets(:)), ...
        double(o.min_length), logical(o.record_parent_ids));
end
