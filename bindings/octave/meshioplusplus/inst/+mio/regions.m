function out = regions(mesh)
%REGIONS  The named regions, an n x 1 struct array with fields name, kind
%   (point, cell or side), dim, tag and entries (num_entries x
%   stride; 1-based, except the facet-ordinal column of a side region).

    out = meshioplusplus_mex('mio_regions', mesh);
end
