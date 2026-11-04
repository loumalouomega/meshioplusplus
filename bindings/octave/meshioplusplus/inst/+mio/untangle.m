function r = untangle(mesh, varargin)
%UNTANGLE  Untangle inverted tetrahedra and hexahedra, keeping connectivity.
%   Options: 'fix_boundary' (true), 'rings' (2), 'smooth' (false) and
%   'max_sweeps' (100). Returns a struct: mesh, num_inverted_before,
%   num_inverted_after, min_quality_before, min_quality_after and sweeps.

    o = mio.internal.opts(varargin, {'fix_boundary', true, 'rings', 2, 'smooth', false, ...
        'max_sweeps', 100}, 'mio.untangle');
    r = meshioplusplus_mex('mio_untangle', mesh, logical(o.fix_boundary), double(o.rings), ...
        logical(o.smooth), double(o.max_sweeps));
end
