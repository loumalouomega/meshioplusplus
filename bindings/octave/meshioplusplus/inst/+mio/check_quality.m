function r = check_quality(mesh, varargin)
%CHECK_QUALITY  Quality gate: test thresholds and inverted/degenerate counts.
%   R = mio.check_quality(MESH, REQUIRE, ...) where REQUIRE is a text of
%   clauses such as 'scaled_jacobian >= 0.2; aspect_ratio <= 5 @ 1%' (''
%   checks only the counts). Options: 'max_inverted' (0) and
%   'max_degenerate' (0); negative disables a count. Returns a struct:
%   passed, num_checks, num_failed, num_cells, num_inverted, num_degenerate
%   and summary (the text both CLIs print).

    o = mio.internal.opts(varargin, {'require', '', 'max_inverted', 0, 'max_degenerate', 0}, ...
        'mio.check_quality');
    r = meshioplusplus_mex('mio_check_quality', mesh, mio.internal.str(o.require, 'require'), ...
        double(o.max_inverted), double(o.max_degenerate));
end
