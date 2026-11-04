function r = project_to_cad(mesh, info, varargin)
%PROJECT_TO_CAD  Move a mesh's points back onto its CAD model.
%   R = mio.project_to_cad(MESH, INFO, ...) moves every point of MESH onto
%   the CAD face or edge its cells name ('cad:face' on triangles and quads,
%   'cad:edge' on lines, as mio.mesh_cad_surface writes them) of the B-rep
%   held by the mio.FormatInfo INFO (roadmap 6.10.1): a point on one edge's
%   lines onto that edge, on two or more (a corner) onto the lowest, a point
%   whose cells name one face onto its surface (Newton point inversion); a
%   point on two faces and no edge stays. Data and regions pass through.
%   Options (positional in this order, or name/value):
%     'face_array'    cell data naming each cell's face ('cad:face')
%     'edge_array'    cell data naming each cell's edge ('cad:edge')
%     'max_distance'  a point farther than this stays; 0 is no limit (0)
%     'trim_samples'  samples per knot span when testing trims (16)
%   R is a struct: mesh, num_projected, num_on_edges, num_on_curves (of
%   those on edges, onto the edge's own 3-D curve), num_corners,
%   num_ambiguous, num_too_far, num_failed, num_outside_trim, max_move and
%   outside_trim (1-based points projected outside their face's trims, kept,
%   never clamped; a column).

    o = mio.internal.opts(varargin, {'face_array', 'cad:face', 'edge_array', 'cad:edge', ...
        'max_distance', 0, 'trim_samples', 16}, 'mio.project_to_cad');
    r = meshioplusplus_mex('mio_project_to_cad', mesh, info, ...
        mio.internal.str(o.face_array, 'face_array'), mio.internal.str(o.edge_array, 'edge_array'), ...
        double(o.max_distance), double(o.trim_samples));
end
