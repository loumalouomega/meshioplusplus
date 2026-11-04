function r = discretize_levelset(mesh, level_set, varargin)
%DISCRETIZE_LEVELSET  Embed a level set conformally with MMG.
%   R = mio.discretize_levelset(MESH, LEVEL_SET, ...) makes the 'isovalue' of
%   the one-component point data LEVEL_SET a conformal interface and remeshes
%   both sides (roadmap 3.2.2, doc/mmg.md). The elements carry cell data
%   'medit:ref' (3 inside, 2 outside) and the cell regions 'inside' and
%   'outside'; the interface is the side region 'interface'. Options:
%     'isovalue'                the level to embed (0)
%     'remove_small_components' MMG's rmc volume fraction, negative to keep
%                               every component (-1)
%     'materials'               [] or a k x 4 matrix with columns ref, split
%                               (0/1), ref_inside and ref_outside: MMG's
%                               multi-material mode ([])
%     'hmin', 'hmax', 'hausd', 'hgrad', 'angle_detection', 'frozen',
%     'required'                as in mio.adapt
%   R is the struct mio.adapt returns.

    o = mio.internal.opts(varargin, {'isovalue', 0, 'remove_small_components', -1, ...
        'materials', [], 'hmin', -1, 'hmax', -1, 'hausd', -1, 'hgrad', -1, ...
        'angle_detection', 45, 'frozen', [], 'required', {}}, 'mio.discretize_levelset');
    materials = o.materials;
    if ~isempty(materials)
        materials = reshape(double(materials), [], 4);
    end
    r = meshioplusplus_mex('mio_discretize_levelset', mesh, ...
        mio.internal.str(level_set, 'level_set'), double(o.isovalue), ...
        double(o.remove_small_components), materials, double(o.hmin), double(o.hmax), ...
        double(o.hausd), double(o.hgrad), double(o.angle_detection), double(o.frozen), ...
        mio.internal.cellstr_arg(o.required, 'required'));
end
