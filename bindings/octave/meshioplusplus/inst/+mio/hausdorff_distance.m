function r = hausdorff_distance(a, b, varargin)
%HAUSDORFF_DISTANCE  Hausdorff distance between the surfaces A and B.
%   Options: 'face_samples' (0: vertices only; k also samples the k^2
%   sub-triangle centroids of every triangle), 'region_a' and 'region_b'
%   ('' for every surface cell). Returns a struct: distance, a_to_b, b_to_a,
%   mean_a_to_b, rms_a_to_b, mean_b_to_a, rms_b_to_a, num_samples_a,
%   num_samples_b, worst_point_a and worst_point_b (1 x 3 rows).

    o = mio.internal.opts(varargin, {'face_samples', 0, 'region_a', '', 'region_b', ''}, ...
        'mio.hausdorff_distance');
    r = meshioplusplus_mex('mio_hausdorff_distance', a, b, double(o.face_samples), ...
        mio.internal.str(o.region_a, 'region_a'), mio.internal.str(o.region_b, 'region_b'));
end
