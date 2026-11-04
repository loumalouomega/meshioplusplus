function n = nurbs_info(info)
%NURBS_INFO  Copy an MFEM NURBS side channel.
%   N = mio.nurbs_info(INFO) copies the patches an MFEM NURBS mesh was
%   sampled from (v16.40.0). N is a struct with
%     patches              n x 1 struct array: degrees (column), knots (a
%                          cell, one column per parametric direction),
%                          control_points (points x space_dimension,
%                          direction 1 fastest), weights (column), fields
%                          (a cell, one points x components matrix of
%                          control values per field)
%     field_names          cellstr
%     topology             struct: dimension, num_vertices, patch_vertices,
%                          patch_attributes, boundary_vertices,
%                          boundary_attributes, edges (one row per edge:
%                          knot vector, vertex, vertex); [] when the source
%                          states none
%     cell_patch, cell_span (cells x 3), boundary_cell_patch
%   Vertices, patches and knot spans are 1-based; attributes and knot-vector
%   ids stay as the file spelled them. The copy remains valid after the
%   handle is released; mio.write_with_info writes the patches back from the
%   handle itself.

    n = meshioplusplus_mex('mio_nurbs_info', info);
end
