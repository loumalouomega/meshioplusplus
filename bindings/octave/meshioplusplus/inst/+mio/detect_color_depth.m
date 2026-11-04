function depth = detect_color_depth(varargin)
%DETECT_COLOR_DEPTH  The colour depth a terminal advertises.
%   DEPTH = mio.detect_color_depth(NO_COLOR, COLOR_TERM, TERM) takes the
%   values of the NO_COLOR, COLORTERM and TERM environment variables, each
%   [] when unset (the default), positionally or as 'no_color',
%   'color_term' and 'term' pairs. DEPTH is 'truecolor', '256', '16' or
%   'mono'. For the current terminal, pass each getenv(...) value, or []
%   where getenv returns '' (unset).

    o = mio.internal.opts(varargin, {'no_color', [], 'color_term', [], 'term', []}, ...
        'mio.detect_color_depth');
    depth = meshioplusplus_mex('mio_detect_color_depth', ...
        mio.internal.str(o.no_color, 'no_color'), ...
        mio.internal.str(o.color_term, 'color_term'), mio.internal.str(o.term, 'term'));
end
