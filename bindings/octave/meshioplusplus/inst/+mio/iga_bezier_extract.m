function m = iga_bezier_extract(info, varargin)
%IGA_BEZIER_EXTRACT  Bezier extraction of a NURBS patch.
%   M = mio.iga_bezier_extract(INFO) extracts patch 1 of the NURBS model in the
%   mio.FormatInfo INFO as VTK Bezier cells; M = mio.iga_bezier_extract(INFO,
%   PATCH) extracts patch PATCH (1-based). Every interior knot is raised to
%   multiplicity p and the net cut at its knot spans: one Bezier curve,
%   quadrilateral or hexahedron cell per span product, (p+1)^d unshared
%   points each, with the cell data HigherOrderDegrees, iga:patch and iga:span
%   and, for a rational patch, the point data RationalWeights. A patch with
%   unequal degrees is refused (raise them with mio.iga_refine first).
%   M is a mio.Mesh.

    if isempty(varargin)
        patch = 1;
    else
        patch = varargin{1};
    end
    m = meshioplusplus_mex('mio_iga_bezier_extract', info, double(patch));
end
