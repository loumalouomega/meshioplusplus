function r = align_datums(mesh, primary, secondary, varargin)
%ALIGN_DATUMS  The 3-2-1 datum alignment of a measured mesh.
%   R = mio.align_datums(MESH, PRIMARY, SECONDARY, ...) moves MESH so the
%   element PRIMARY sets Z and SECONDARY X (roadmap 7.3.4). Option 'origin'
%   ('': the primary's point) names a point-like element or a third plane. R
%   is a struct: mesh and the 4 x 4 matrix.

    o = mio.internal.opts(varargin, {'origin', ''}, 'mio.align_datums');
    r = meshioplusplus_mex('mio_align_datums', mesh, mio.internal.str(primary, 'primary'), ...
        mio.internal.str(secondary, 'secondary'), o.origin);
end
