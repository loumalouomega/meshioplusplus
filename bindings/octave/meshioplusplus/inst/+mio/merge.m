function m = merge(meshes, varargin)
%MERGE  Concatenate several meshes into one new mesh.
%   M = mio.merge(MESHES) where MESHES is a non-empty cell array of
%   mio.MeshHandle objects. Options (positional or name/value):
%     'weld'                 merge coincident points (false)
%     'atol'                 weld tolerance (1e-12)
%     'source_tag'           add an integer 'source_mesh_id' cell-data array (true)
%     'data_policy'          'intersection' (keep keys present in every input)
%                            or 'fill' (union, NaN for missing rows)
%     'drop_duplicate_cells' drop cells repeated across inputs (false)

    o = mio.internal.opts(varargin, {'weld', false, 'atol', 1e-12, 'source_tag', true, ...
        'data_policy', 'intersection', 'drop_duplicate_cells', false}, 'mio.merge');
    if ~iscell(meshes)
        meshes = num2cell(meshes);
    end
    policy = mio.internal.choice(o.data_policy, {'intersection', 'fill'}, 'data_policy');
    m = meshioplusplus_mex('mio_merge', meshes, logical(o.weld), double(o.atol), ...
        logical(o.source_tag), policy, logical(o.drop_duplicate_cells));
end
