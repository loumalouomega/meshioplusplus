function provenance_target(format, varargin)
%PROVENANCE_TARGET  Record what was actually written, on the active provenance scope.
%   mio.provenance_target(FORMAT, ...) records the target format and the
%   encoding / codec / float format actually used, as opposed to what was
%   requested. Options 'encoding', 'codec' and 'float_format' ('' when not
%   applicable, the default). A no-op outside a scope (see
%   mio.with_provenance).

    o = mio.internal.opts(varargin, {'encoding', '', 'codec', '', 'float_format', ''}, ...
        'mio.provenance_target');
    meshioplusplus_mex('mio_provenance_set_target', mio.internal.str(format, 'format'), ...
        mio.internal.str(o.encoding, 'encoding'), mio.internal.str(o.codec, 'codec'), ...
        mio.internal.str(o.float_format, 'float_format'));
end
