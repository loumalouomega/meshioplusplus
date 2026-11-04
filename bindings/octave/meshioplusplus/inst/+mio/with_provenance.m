function varargout = with_provenance(fn, varargin)
%WITH_PROVENANCE  Call a function with a provenance scope open.
%   [...] = mio.with_provenance(FN) calls the function handle FN (taking no
%   arguments) with a provenance scope open and returns FN's outputs. Any
%   mio.write inside FN records the richer provenance block -- source,
%   target, operation chain, conversion assumptions, timestamp -- into the
%   file's header slot, where the format has one. With no scope open a
%   writer emits only the one-line credit, which is the default everywhere.
%
%   This takes a function handle where the R binding's mio_with_provenance
%   takes an unevaluated expression. The scope is closed on exit even if FN
%   throws (an onCleanup object pairs the underlying begin/end C calls).
%
%   Option 'mode': 'best_effort' (default), 'off' or 'required'. 'required'
%   errors for a format whose header slot cannot hold a block at all; it
%   does not error merely because a slot is too small for the full block,
%   since degrading to the credit line is that slot's honest maximum. See
%   doc/provenance.md.
%
%   Example:
%       mio.with_provenance(@() write_it(mesh));
%   where write_it calls
%       mio.provenance_source('bracket.msh', 'gmsh');
%       mio.provenance_note('regions-dropped', 'Side regions have no VTU equivalent');
%       mio.write(mesh, 'out.vtu');

    o = mio.internal.opts(varargin, {'mode', 'best_effort'}, 'mio.with_provenance');
    if ~isa(fn, 'function_handle')
        error('meshioplusplus:argument', '`fn` must be a function handle');
    end
    mode = mio.internal.choice(o.mode, {'off', 'best_effort', 'required'}, 'mode');
    meshioplusplus_mex('mio_provenance_begin', mode);
    cleanup = onCleanup(@() meshioplusplus_mex('mio_provenance_end'));
    if nargout > 0
        [varargout{1:nargout}] = fn();
    else
        fn();
    end
    clear cleanup;
end
