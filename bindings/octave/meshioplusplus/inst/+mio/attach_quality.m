function m = attach_quality(mesh)
%ATTACH_QUALITY  A copy of MESH carrying the cell-quality metrics.
%   Q = mio.attach_quality(MESH) returns a new mesh with one 'quality:<metric>'
%   cell-data array per quality metric.

    m = meshioplusplus_mex('mio_attach_quality', mesh);
end
