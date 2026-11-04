function sequence_to_timeseries(seq, out_path, varargin)
%SEQUENCE_TO_TIMESERIES  Fan-in: write every step into one multi-step file.
%   mio.sequence_to_timeseries(SEQ, OUT_PATH, ...) streams the steps of SEQ
%   -- one mesh alive at a time -- into OUT_PATH. A format that cannot hold a
%   series fails naming itself and pointing at {step}, never a silent
%   truncation to the first step.
%   Options (positional in this order, or name/value):
%     'out_format'  forced output format, [] to infer
%     'ascii'       select XDMF's "XML" data format (everything inline, no
%                   HDF5 needed) instead of the default "HDF" (false)

    o = mio.internal.opts(varargin, {'out_format', [], 'ascii', false}, ...
        'mio.sequence_to_timeseries');
    meshioplusplus_mex('mio_sequence_to_timeseries', seq, mio.internal.str(out_path, 'out_path'), ...
        mio.internal.str(o.out_format, 'out_format'), logical(o.ascii));
end
