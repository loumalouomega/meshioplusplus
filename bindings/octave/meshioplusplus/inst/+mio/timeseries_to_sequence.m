function timeseries_to_sequence(in_path, out_pattern, varargin)
%TIMESERIES_TO_SEQUENCE  Fan-out: write each step of a multi-step file.
%   mio.timeseries_to_sequence(IN_PATH, OUT_PATTERN, ...) writes every step
%   of the multi-step file IN_PATH to OUT_PATTERN, which must contain {step}
%   or {index}.
%   Options (positional in this order, or name/value):
%     'in_format'   forced input format, [] to infer
%     'out_format'  forced output format, [] to infer

    o = mio.internal.opts(varargin, {'in_format', [], 'out_format', []}, ...
        'mio.timeseries_to_sequence');
    meshioplusplus_mex('mio_timeseries_to_sequence', mio.internal.str(in_path, 'in_path'), ...
        mio.internal.str(o.in_format, 'in_format'), mio.internal.str(out_pattern, 'out_pattern'), ...
        mio.internal.str(o.out_format, 'out_format'));
end
