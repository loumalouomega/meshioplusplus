function r = stats(mesh)
%STATS  Summary statistics of a mesh.
%   R = mio.stats(MESH) is a struct: num_points, num_cells, bbox_min,
%   bbox_max, extent and centroid (each 1 x 3), total_area, signed_volume,
%   unsigned_volume and num_inverted.

    r = meshioplusplus_mex('mio_stats', mesh);
end
