function mesh = add_region(mesh, name, kind, entries, varargin)
%ADD_REGION  Add a named region.
%   KIND is 'point', 'cell' or 'side'. ENTRIES is a vector of 1-based point
%   or global cell indices, or for a side region an n x 2 matrix of [cell,
%   facet] rows (the facet ordinal is NOT shifted). Options: 'dim' (-1),
%   'tag' (-1).

    o = mio.internal.opts(varargin, {'dim', -1, 'tag', -1}, 'mio.add_region');
    meshioplusplus_mex('mio_add_region', mesh, mio.internal.str(name, 'name'), ...
        mio.internal.str(kind, 'kind'), double(entries), double(o.dim), double(o.tag));
end
