function r = constrain_primitives(mesh, varargin)
%CONSTRAIN_PRIMITIVES  Fit related primitives of a fitted skin together.
%   R = mio.constrain_primitives(MESH, ...) relates coaxial, parallel,
%   perpendicular and equal-radius primitives (roadmap 7.7.3). Options:
%   'face_tolerance' (0), 'angle_tolerance', 'distance_tolerance' and
%   'min_fraction' (-1: the pass's own), and 'force' and 'forbid' as n x 3
%   cell arrays {patch_a, patch_b, relation} (relation 'coaxial', 'parallel',
%   'perpendicular' or 'equal_radius'). R is a struct: mesh,
%   num_constrained, num_groups, face_tolerance, and per relation applied
%   patch_a, patch_b, relation, angle_deg, distance and forced; unmet holds
%   the forced relations not kept, as rows {patch_a, patch_b, relation}.

    o = mio.internal.opts(varargin, {'face_tolerance', 0, 'angle_tolerance', -1, ...
        'distance_tolerance', -1, 'min_fraction', -1, 'force', {}, 'forbid', {}}, ...
        'mio.constrain_primitives');
    kinds = {'coaxial', 'parallel', 'perpendicular', 'equal_radius'};
    r = meshioplusplus_mex('mio_constrain_primitives', mesh, double(o.face_tolerance), ...
        double(o.angle_tolerance), double(o.distance_tolerance), double(o.min_fraction), ...
        triples(o.force, kinds), triples(o.forbid, kinds));
    r.relation = kinds(r.relation + 1);
    r.relation = r.relation(:);
    r.forced = logical(r.forced);
    u = r.unmet;
    r.unmet = cell(size(u, 1), 3);
    for i = 1:size(u, 1)
        r.unmet(i, :) = {u(i, 1), u(i, 2), kinds{u(i, 3) + 1}};
    end
end

function t = triples(list, kinds)
    t = zeros(1, 3 * size(list, 1));
    for i = 1:size(list, 1)
        k = mio.internal.choice(list{i, 3}, kinds, 'relation');
        t(3 * i - 2:3 * i) = [double(list{i, 1}), double(list{i, 2}), k];
    end
end
