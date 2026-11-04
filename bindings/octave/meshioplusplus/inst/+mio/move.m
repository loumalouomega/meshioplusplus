function r = move(mesh, displacement, moving, varargin)
%MOVE  Move a mesh by a boundary displacement with MMG.
%   R = mio.move(MESH, DISPLACEMENT, MOVING, ...): the boundary facets of the
%   regions MOVING (a cellstr of cell or side region names) follow the
%   point-data vector DISPLACEMENT, MMG extends it inward by linear
%   elasticity (ELAS) and remeshes as the mesh deforms; the rest of the
%   boundary stays in place (roadmap 3.2.4, doc/mmg.md). Fails naming
%   MESHIOPLUSPLUS_WITH_MMG without MMG, or ELAS without it
%   (mio.has_mmg_elas). Options:
%     'lagrangian' MMG's lag: 0 moves only, 1 also swaps, 2 also inserts (1)
%     'hmin', 'hmax', 'hausd', 'hgrad', 'angle_detection', 'frozen',
%     'required'   as in mio.adapt
%   R is the struct mio.adapt returns.

    o = mio.internal.opts(varargin, {'lagrangian', 1, 'hmin', -1, 'hmax', -1, 'hausd', -1, ...
        'hgrad', -1, 'angle_detection', 45, 'frozen', [], 'required', {}}, 'mio.move');
    r = meshioplusplus_mex('mio_move', mesh, mio.internal.str(displacement, 'displacement'), ...
        mio.internal.cellstr_arg(moving, 'moving'), double(o.lagrangian), double(o.hmin), ...
        double(o.hmax), double(o.hausd), double(o.hgrad), double(o.angle_detection), ...
        double(o.frozen), mio.internal.cellstr_arg(o.required, 'required'));
end
