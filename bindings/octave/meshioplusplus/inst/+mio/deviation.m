function r = deviation(scan, varargin)
%DEVIATION  Signed deviation of a scan from a reference.
%   R = mio.deviation(SCAN, ...) measures the signed distance of every point
%   of SCAN to the surface 'reference' (positive outside it, excess material),
%   or, with 'reference' [], to the fitted primitive 'type' with its ten
%   fit:patch_params entries 'params' (roadmap 7.4). Options:
%     'reference'        the reference surface, or [] for a primitive ([])
%     'type'             'plane', 'sphere', 'cylinder', 'cone' or 'torus' ('')
%     'params'           the primitive's ten parameters ([])
%     'tolerance'        report the share within +-tolerance, 0 for none (0)
%     'bins'             histogram bins (20)
%     'range'            [lo hi] of the histogram, [] for +-max|d| ([])
%     'max_normal_angle' the far-side rule's angle in degrees, 180 off (90)
%     'reference_region' named cell region of the reference ('')
%     'register'         register the scan to the reference first (false);
%                        mesh is then the registered scan
%   R is a struct: mesh (with the 'deviation:distance' point data),
%   num_points, num_measured, min, max, mean, rms, sigma, percentiles
%   (1 x 5: p1, p5, p50, p95, p99), histogram_edges, histogram_counts,
%   num_outside_range, share_within, suggested_range, num_without_normal,
%   warnings (n x 1 cellstr), registered, matrix (4x4) and registration_rms.

    o = mio.internal.opts(varargin, {'reference', [], 'type', '', 'params', [], ...
        'tolerance', 0, 'bins', 20, 'range', [], 'max_normal_angle', 90, ...
        'reference_region', '', 'register', false}, 'mio.deviation');
    types = {'plane', 'sphere', 'cylinder', 'cone', 'torus'};
    ptype = 0;
    if isempty(o.reference)
        ptype = find(strcmp(types, o.type));
        if isempty(ptype)
            error('mio:deviation', 'mio.deviation: `type` must name a primitive');
        end
    end
    rng = [0 0];
    if ~isempty(o.range)
        rng = double(o.range(:)');
    end
    r = meshioplusplus_mex('mio_deviation', scan, o.reference, double(ptype), ...
        double(o.params), double(o.tolerance), double(o.bins), rng(1), rng(2), ...
        double(o.max_normal_angle), mio.internal.str(o.reference_region, 'reference_region'), ...
        logical(o.register));
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
