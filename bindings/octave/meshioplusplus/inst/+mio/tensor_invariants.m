function out = tensor_invariants(mesh, location, varargin)
%TENSOR_INVARIANTS  Von Mises, principal, hydrostatic and deviatoric invariants.
%   OUT = mio.tensor_invariants(MESH, LOCATION) adds the invariants of the
%   symmetric (6-component, xx yy zz xy yz zx) or general 3x3 (9-component,
%   row-major) tensor arrays at LOCATION ('point' or 'cell'; 'field' is
%   refused) to a new mesh, named prefix + <name> + '_mises' + suffix, and
%   likewise _principal (eigenvalues, ascending), _hydrostatic and
%   _deviatoric. Options:
%     'names'      tensor arrays ([] or {} for every 6- or 9-component array)
%     'outputs'    any of 'mises', 'principal', 'hydrostatic', 'deviatoric'
%                  ([] or {} for all four)
%     'prefix'     prepended to every output name ([] for none)
%     'suffix'     appended after the invariant's name ([] for none)
%     'overwrite'  replace existing output arrays (true)

    o = mio.internal.opts(varargin, {'names', [], 'outputs', [], 'prefix', [], ...
        'suffix', [], 'overwrite', true}, 'mio.tensor_invariants');
    flags = {'mises', 'principal', 'hydrostatic', 'deviatoric'};
    outputs = mio.internal.cellstr_arg(o.outputs, 'outputs');
    if isempty(outputs)
        mask = 15;
    else
        mask = 0;
        for k = 1:numel(outputs)
            hit = find(strcmp(outputs{k}, flags), 1);
            if isempty(hit)
                error('meshioplusplus:argument', ['`outputs` entries must be ''mises'', ' ...
                      '''principal'', ''hydrostatic'' or ''deviatoric'', got ''%s'''], outputs{k});
            end
            mask = bitor(mask, 2 ^ (hit - 1));
        end
    end
    out = meshioplusplus_mex('mio_tensor_invariants', mesh, ...
        mio.internal.choice(location, {'point', 'cell', 'field'}, 'location'), ...
        mio.internal.cellstr_arg(o.names, 'names'), double(mask), ...
        mio.internal.str(o.prefix, 'prefix'), mio.internal.str(o.suffix, 'suffix'), ...
        logical(o.overwrite));
end
