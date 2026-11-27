function r = streamlines(mesh, array, varargin)
%STREAMLINES  Integrate streamlines of a vector point-data field.
%   Options: 'seeds' (40), 'direction' ('both'|'forward'|'backward'),
%   'integrator' ('rk4'|'rk45'), 'max_length' (0.5), 'step' (0.0),
%   'max_steps' (0), 'tolerance' (1e-6).

    o = mio.internal.opts(varargin, {'seeds', 40, 'direction', 'both', ...
        'integrator', 'rk4', 'max_length', 0.5, 'step', 0.0, ...
        'max_steps', 0, 'tolerance', 1e-6}, 'mio.streamlines');
    switch o.direction
        case 'both'
            d = 0;
        case 'forward'
            d = 1;
        case 'backward'
            d = 2;
        otherwise
            error('mio:streamlines', 'unknown direction ''%s''', o.direction);
    end
    switch o.integrator
        case 'rk4'
            ig = 0;
        case 'rk45'
            ig = 1;
        otherwise
            error('mio:streamlines', 'unknown integrator ''%s''', o.integrator);
    end
    r = meshioplusplus_mex('mio_streamlines', mesh, char(array), ...
        double(o.seeds), double(d), double(ig), double(o.max_length), ...
        double(o.step), double(o.max_steps), double(o.tolerance));
end
