function m = med_read_named(path, name, varargin)
%MED_READ_NAMED  Read the MED mesh called NAME. Options: 'time_step' (0), 'lenient' (false).

    o = mio.internal.opts(varargin, {'time_step', 0, 'lenient', false}, 'mio.med_read_named');
    m = meshioplusplus_mex('mio_med_read_named', mio.internal.str(path, 'path'), ...
        mio.internal.str(name, 'name'), double(o.time_step), logical(o.lenient));
end
