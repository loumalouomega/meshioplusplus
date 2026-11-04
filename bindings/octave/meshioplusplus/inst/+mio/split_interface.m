function r = split_interface(mesh, side, varargin)
%SPLIT_INTERFACE  Split MESH along the named Side region SIDE.
%   Duplicates the incident-cell fans along SIDE, preserving data and
%   regions. Options: 'add_cohesive' (false: insert cohesive cells on
%   two-owner facets). Returns a struct: mesh, num_duplicated_points and
%   num_cohesive_cells.

    o = mio.internal.opts(varargin, {'add_cohesive', false}, 'mio.split_interface');
    r = meshioplusplus_mex('mio_split_interface', mesh, mio.internal.str(side, 'side'), ...
        logical(o.add_cohesive));
end
