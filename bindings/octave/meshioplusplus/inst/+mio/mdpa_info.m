function d = mdpa_info(info)
%MDPA_INFO  Copy an MDPA (Kratos) side channel into plain values.
%   D = mio.mdpa_info(INFO) copies the side channel of a mio.FormatInfo
%   read from an .mdpa file. D is a struct with fields
%     properties       n x 1 struct array: id, values (a struct)
%     entity_names     n x 1 struct array: name, is_condition
%     skipped          n x 1 cellstr
%     model_part_data  struct of the non-numeric entries (key -> value)
%     tables           n x 1 struct array: header, values (rows x columns)
%     geometries       n x 1 struct array: name, type, connectivity
%                      (num_geometries x nodes_per_geometry, 1-based), ids
%     mesh_blocks      n x 1 struct array: id, data, nodes (1-based),
%                      element_ids, condition_ids
%     submodelparts    n x 1 struct array: name, data, tables,
%                      geometry_ids, constraint_ids
%     raw_blocks       n x 1 struct array: header, body, terminator
%   Id vectors are columns; file ids stay as the file spelled them (as
%   double). Data entries become struct fields named by their key (a key that
%   is not a valid identifier is made one, a repeated key gets a _2 suffix).
%   The copy remains valid after the handle is released.

    d = meshioplusplus_mex('mio_mdpa_info', info);
end
