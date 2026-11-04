function m = crop_bbox(mesh, lo, hi, varargin)
%CROP_BBOX  Keep the cells inside the box [LO, HI]; returns a new mesh.
%   LO and HI are 3-vectors. Options: 'mode' ('all': every node inside, or
%   'any': at least one) and 'record_ids' (false: attach the 'crop:'
%   original-id arrays).

    o = mio.internal.opts(varargin, {'mode', 'all', 'record_ids', false}, 'mio.crop_bbox');
    m = meshioplusplus_mex('mio_crop_bbox', mesh, double(lo), double(hi), ...
        mio.internal.choice(o.mode, {'all', 'any'}, 'mode'), logical(o.record_ids));
end
