function provenance_note(category, detail)
%PROVENANCE_NOTE  Record one conversion assumption against the active provenance scope.
%   mio.provenance_note(CATEGORY, DETAIL) -- CATEGORY a short stable slug,
%   e.g. 'regions-dropped' or 'dtype'; DETAIL human-readable specifics. A
%   no-op outside a scope (see mio.with_provenance), so it is safe to call
%   unconditionally.

    meshioplusplus_mex('mio_provenance_note', mio.internal.str(category, 'category'), ...
        mio.internal.str(detail, 'detail'));
end
