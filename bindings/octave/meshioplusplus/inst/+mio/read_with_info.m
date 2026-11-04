function r = read_with_info(path, varargin)
%READ_WITH_INFO  Read a mesh and keep the format's side channel.
%   R = mio.read_with_info(PATH, ...) reads PATH like mio.read and keeps
%   what the format carries that a mesh cannot hold -- MDPA tables,
%   geometries, Mesh blocks, sub-model-part data and blocks such as
%   Constraints; Gmsh periodic links and bounding entities; an MFEM NURBS
%   model; a CAD JSON / STEP / IGES B-rep -- in a mio.FormatInfo handle that
%   mio.write_with_info puts back. mio.read refuses such a file by name (or,
%   with 'lenient', skips the content); Gmsh periodic files need this
%   channel even with 'lenient'.
%
%   Options (positional in this order, or name/value):
%     'format'              explicit format name ([] to infer)
%     'lenient'             skip (with a warning) what not even the side
%                           channel holds (false)
%     'deflection'          STEP/IGES: chordal deflection of the
%                           triangulation, a length; <= 0 is automatic (0)
%     'angular_deflection'  STEP/IGES: angular deflection in radians (0.5)
%     'heal'                STEP/IGES: heal the shape before it is extracted
%                           (ShapeFix_Shape, then sewing); the report is the
%                           "heal" member of mio.brep_info(info).extra (false)
%     'sew_tolerance'       STEP/IGES: sewing tolerance; <= 0 is automatic
%                           when healing, a positive value sews without
%                           'heal' (0)
%
%   R is a struct with fields 'mesh' (a mio.MeshHandle) and 'info' (a
%   mio.FormatInfo, or [] for a format with no side channel). The handle is
%   freed when cleared, or now by mio.format_info_release.

    o = mio.internal.opts(varargin, {'format', [], 'lenient', false, 'deflection', 0, ...
        'angular_deflection', 0.5, 'heal', false, 'sew_tolerance', 0}, 'mio.read_with_info');
    r = meshioplusplus_mex('mio_read_with_info', mio.internal.str(path, 'path'), ...
        mio.internal.str(o.format, 'format'), logical(o.lenient), double(o.deflection), ...
        double(o.angular_deflection), logical(o.heal), double(o.sew_tolerance));
end
