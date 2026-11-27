function r = transfer_operator(coarse, fine, varargin)
%TRANSFER_OPERATOR  Sparse CSR prolongation/restriction operator between meshes.
%   Options: 'mode' ('auto'|'nested'|'barycentric'), 'transpose' (false),
%   'tolerance' (1e-12), 'fallback_nearest' (false).

    o = mio.internal.opts(varargin, {'mode', 'auto', 'transpose', false, ...
        'tolerance', 1e-12, 'fallback_nearest', false}, 'mio.transfer_operator');
    switch o.mode
        case 'auto'
            m = 0;
        case 'nested'
            m = 1;
        case 'barycentric'
            m = 2;
        otherwise
            error('mio:transfer_operator', 'unknown mode ''%s''', o.mode);
    end
    r = meshioplusplus_mex('mio_transfer_operator', coarse, fine, double(m), ...
        logical(o.transpose), double(o.tolerance), logical(o.fallback_nearest));
end
