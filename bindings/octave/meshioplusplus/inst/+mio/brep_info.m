function b = brep_info(info)
%BREP_INFO  Copy a B-rep side channel.
%   B = mio.brep_info(INFO) copies the trimmed NURBS faces, loops, edges and
%   vertices of a Kratos CAD JSON handle (or a STEP/IGES one read through
%   OpenCASCADE). B is a struct with
%     patches   struct array: degrees (2 x 1), knots (2 x 1 cell, full
%               form), control_points (points x 3, direction 1 fastest),
%               weights
%     breps     struct array: id, name, faces, edges, vertices (1-based),
%               extra
%     faces     struct array: id, name, patch (1-based), swapped_normal,
%               trimmed, rational, reduced_knots, control_point_ids, loops
%               (struct array: type, trims, extra), embedded_edges,
%               surface_extra, extra; each trim: trim_index, direction,
%               curve, extra
%     edges     struct array: id, name, curve (or []), topology, extra
%     vertices  struct array: id, name, topology, extra
%     extra     the file's other members
%   Each topology entry has id, name, trim_index, relative_direction,
%   local_coordinates, extra. A curve has degree, knots, reduced_knots,
%   rational, control_points, weights, ids, active_range, extra. Ids are -1
%   where the file names the entity instead; the extra members are compact
%   JSON text. The copy remains valid after the handle is released;
%   mio.write_with_info writes the B-rep back as CAD JSON from the handle.

    b = meshioplusplus_mex('mio_brep_info', info);
end
