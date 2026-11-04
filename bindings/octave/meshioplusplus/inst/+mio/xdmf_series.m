function s = xdmf_series(path, varargin)
%XDMF_SERIES  Open a transient (time-series) XDMF writer.
%   S = mio.xdmf_series(PATH) starts a series at PATH (.xdmf / .xmf) and
%   returns a mio.XdmfSeries. It is the one writer mio.write cannot express:
%   the grid is written once with mio.xdmf_series_write_points_cells and each
%   solve appends a cheap step with mio.xdmf_series_write_data.
%
%   The .xdmf light data is buffered until the series is flushed or
%   finalized, so the file is only readable after mio.xdmf_series_flush,
%   mio.xdmf_series_finalize or mio.xdmf_series_release. Heavy data for
%   'HDF' goes to a <path minus extension>.h5 sibling of the .xdmf.
%
%   The handle is freed (and the file finalized) when its last reference
%   goes away; mio.xdmf_series_release does it immediately and is
%   idempotent. A write failure during that implicit finalize cannot be
%   reported, which is why mio.xdmf_series_finalize exists.
%
%   Options, positional in this order or as name/value pairs:
%     'data_format'  'HDF' (default; needs an HDF5-enabled library), 'XML'
%                    (everything inline in the .xdmf) or 'Binary'. An
%                    unknown format, or 'HDF' without HDF5, is an error.
%     'gzip_level'   gzip level for 'HDF' datasets; negative (-1, the
%                    default) means no compression
%     'mode'         'truncate' (default) starts a fresh series; 'append'
%                    continues the one already at PATH, if any (a path with
%                    no file yet is simply a fresh series)
%     'auto_flush'   flush after every mio.xdmf_series_write_data (false).
%                    Off by default because a flush re-serializes the whole
%                    document, making per-step flushing quadratic.
%
%   Example:
%     s = mio.xdmf_series('simulation.xdmf', 'XML');
%     mio.xdmf_series_write_points_cells(s, mesh);
%     for k = 0:9
%       mio.xdmf_series_write_data(s, k * 0.1, mesh);
%     end
%     mio.xdmf_series_finalize(s);
%     mio.xdmf_series_release(s);

    o = mio.internal.opts(varargin, {'data_format', 'HDF', 'gzip_level', -1, ...
        'mode', 'truncate', 'auto_flush', false}, 'mio.xdmf_series');
    modes = {'truncate', 'append'};
    mode = modes{mio.internal.choice(o.mode, modes, 'mode') + 1};
    s = meshioplusplus_mex('mio_xdmf_series_create', mio.internal.str(path, 'path'), ...
        mio.internal.str(o.data_format, 'data_format'), double(o.gzip_level), mode, ...
        logical(o.auto_flush));
end
