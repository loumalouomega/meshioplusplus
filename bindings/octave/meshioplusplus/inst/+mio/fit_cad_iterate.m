function r = fit_cad_iterate(scan, varargin)
%FIT_CAD_ITERATE  Segment, fit and build a solid in rounds towards a target.
%   R = mio.fit_cad_iterate(SCAN, ...) (roadmap 7.7.2). Options: 'step_path'
%   (''), 'part_class' ('any', 'analytic', 'cast' or 'freeform'),
%   'target_within' and 'target_explained' (0), 'segment_tolerance' and
%   'fit_tolerance' (0: the class's), 'max_rounds' (4), 'worst_faces' (3),
%   'refine_factor' (0.25) and 'backend' ('auto', 'occt' or 'native'). R is a
%   struct: mesh (the best solid), num_rounds, best_round, target_met,
%   stop_reason, target_within, share_within, explained, hausdorff_max,
%   num_faces, round_share_within and round_explained.

    o = mio.internal.opts(varargin, {'step_path', '', 'part_class', 'any', ...
        'target_within', 0, 'target_explained', 0, 'segment_tolerance', 0, ...
        'fit_tolerance', 0, 'max_rounds', 4, 'worst_faces', 3, 'refine_factor', 0.25, ...
        'backend', 'auto'}, 'mio.fit_cad_iterate');
    pc = mio.internal.choice(o.part_class, {'any', 'analytic', 'cast', 'freeform'}, 'part_class');
    be = mio.internal.choice(o.backend, {'auto', 'occt', 'native'}, 'backend');
    r = meshioplusplus_mex('mio_fit_cad_iterate', scan, mio.internal.str(o.step_path, 'step_path'), ...
        double(pc), double(o.target_within), double(o.target_explained), ...
        double(o.segment_tolerance), double(o.fit_tolerance), double(o.max_rounds), ...
        double(o.worst_faces), double(o.refine_factor), double(be));
    r.target_met = logical(r.target_met);
    reasons = {'target', 'no improvement', 'rounds', 'nothing to refine'};
    r.stop_reason = reasons{r.stop_reason + 1};
end
