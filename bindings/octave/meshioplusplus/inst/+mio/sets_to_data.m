function out = sets_to_data(mesh, location, varargin)
%SETS_TO_DATA  Convert point/cell region-backed sets into integer labels.
%   OUT = mio.sets_to_data(MESH, LOCATION) returns a new mesh in which the
%   sets at LOCATION ('point' or 'cell') become one scalar integer data
%   array, and the sets are removed. Labels are zero-based set ordinals (not
%   indices): later overlapping sets win and uncovered entities are -1. Side
%   regions, geometry, unrelated data and property sets are preserved.
%   Options:
%     'data_name'  output array name ([] joins the set names)
%     'join_char'  separator used when joining set names ('-')
%     'order'      cellstr naming every set exactly once ([] for the
%                  native region-name order)
%   See also mio.data_to_sets.

    o = mio.internal.opts(varargin, {'data_name', [], 'join_char', '-', 'order', []}, ...
        'mio.sets_to_data');
    out = meshioplusplus_mex('mio_sets_to_data', mesh, ...
        mio.internal.choice(location, {'point', 'cell', 'field'}, 'location'), ...
        mio.internal.str(o.data_name, 'data_name'), mio.internal.str(o.join_char, 'join_char'), ...
        mio.internal.cellstr_arg(o.order, 'order'));
end
