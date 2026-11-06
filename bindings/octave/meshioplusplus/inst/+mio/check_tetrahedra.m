function r = check_tetrahedra(mesh, varargin)
%CHECK_TETRAHEDRA  Whether the tetrahedra of a mesh are usable.
%   R = mio.check_tetrahedra(MESH, ...) examines the tetra and tetra10 cells (a
%   quadratic tetrahedron by its four corners) and counts the tetrahedra that
%   are inverted (a signed volume below minus the tolerance) or degenerate (a
%   volume within it), and the pairs that overlap across a shared face (no
%   other tetrahedron uses it, and their opposite nodes lie on the same side of
%   it). The mesh is not changed. Options:
%     'zero_volume_factor'  the volume tolerance as a fraction of the longest
%                           edge cubed, plus the effect of six significant
%                           digits of written coordinates (1e-6)
%   R is a struct: num_tetrahedra, num_inverted, num_degenerate, num_overlaps,
%   min_volume, max_volume, valid (logical: no inverted or degenerate
%   tetrahedron and no overlap) and warnings (n x 1 cellstr).

    o = mio.internal.opts(varargin, {'zero_volume_factor', []}, 'mio.check_tetrahedra');
    if isempty(o.zero_volume_factor)
        factor = [];
    else
        factor = double(o.zero_volume_factor);
    end
    r = meshioplusplus_mex('mio_check_tetrahedra', mesh, factor);
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
