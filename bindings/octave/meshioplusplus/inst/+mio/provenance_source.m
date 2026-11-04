function provenance_source(path, format)
%PROVENANCE_SOURCE  Record where a mesh came from, on the active provenance scope.
%   mio.provenance_source(PATH, FORMAT) records the source file path and
%   format name. A no-op outside a scope (see mio.with_provenance).

    meshioplusplus_mex('mio_provenance_set_source', mio.internal.str(path, 'path'), ...
        mio.internal.str(format, 'format'));
end
