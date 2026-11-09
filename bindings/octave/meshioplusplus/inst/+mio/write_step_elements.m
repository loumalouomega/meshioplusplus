function r = write_step_elements(path, mesh, varargin)
%WRITE_STEP_ELEMENTS  Fitted primitives as STEP elements.
%   R = mio.write_step_elements(PATH, MESH, ...) writes the fitted primitives
%   of MESH (mio.segment, then mio.fit_primitives) to the STEP file PATH as
%   exact geometry a CAD user can snap to: bounded planes, balls, lidded
%   cylinders, frusta and tori (roadmap 7.6.5). Options: 'patches' ([]: every
%   fitted patch) and 'segments' (64). R is a struct: num_planes, num_spheres,
%   num_cylinders, num_cones, num_tori, num_skipped and warnings.

    o = mio.internal.opts(varargin, {'patches', [], 'segments', 64}, 'mio.write_step_elements');
    r = meshioplusplus_mex('mio_write_step_elements', path, mesh, double(o.patches(:)), ...
        double(o.segments));
    if isempty(r.warnings)
        r.warnings = cell(0, 1);
    else
        r.warnings = strsplit(r.warnings, char(10));
        r.warnings = r.warnings(:);
    end
end
