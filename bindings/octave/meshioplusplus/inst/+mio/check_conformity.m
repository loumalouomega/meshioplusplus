function r = check_conformity(mesh, varargin)
%CHECK_CONFORMITY  Whether a mesh is conforming.
%   R = mio.check_conformity(MESH, ...) checks the mesh's volume cells (else
%   its 2-D cells) for hanging nodes, T-junctions, and partly overlapping or
%   unwelded pairs of unmatched facets (roadmap 3.2.1). Options:
%     'tolerance'    distances at or below this fraction of the examined
%                    facets' bounding-box diagonal count as zero (1e-9)
%     'attach_flags' attach the 'conformity:hanging' point data (1 hanging
%                    node, 2 T-junction) (true)
%   R is a struct: mesh, dimension, num_hanging_nodes, num_t_junctions,
%   num_partial_overlaps, num_unwelded_facets, num_boundary_facets,
%   num_non_manifold_facets, conforming (logical), interface_mode
%   ('conforming', 'weld' or 'proximity') and warnings (n x 1 cellstr).

    o = mio.internal.opts(varargin, {'tolerance', 1e-9, 'attach_flags', true}, ...
        'mio.check_conformity');
    r = meshioplusplus_mex('mio_check_conformity', mesh, double(o.tolerance), ...
        logical(o.attach_flags));
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
