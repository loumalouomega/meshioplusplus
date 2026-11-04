function c = quality_counts(mesh)
%QUALITY_COUNTS  Count the inverted and degenerate cells of MESH.
%   C = mio.quality_counts(MESH) returns a struct with fields num_cells,
%   num_inverted and num_degenerate.

    c = meshioplusplus_mex('mio_quality_counts', mesh);
end
