function r = mass_properties(mesh, varargin)
%MASS_PROPERTIES  Mass, centre of mass and inertia of a body.
%   Integrates a density over the volume cells of MESH, over the body a closed
%   triangle or polygon surface encloses, or (with 'thickness') over a shell.
%   Options: 'density' (1), 'density_array' ('': the name of a scalar
%   cell_data array of per-cell densities), 'region' ('': restrict the total to
%   this named cell region), 'thickness' (0: above zero, surface cells are a
%   shell of this thickness) and 'about' ('centroid' or 'origin': the point the
%   inertia tensor is taken about).
%   Returns a struct: kind ('solid' or 'shell'), total (a struct of name,
%   num_cells, measure, mass, centre (1 x 3), inertia (1 x 6, xx yy zz xy yz
%   zx), principal_moments (1 x 3, ascending), principal_axes (3 x 3, one axis
%   per column) and gyration_radii (1 x 3)), regions (a struct array of the
%   same fields, one per named cell region) and warnings (a cell array of
%   char).

    o = mio.internal.opts(varargin, {'density', 1, 'density_array', '', 'region', '', ...
        'thickness', 0, 'about', 'centroid'}, 'mio.mass_properties');
    about = mio.internal.str(o.about, 'about');
    if strcmp(about, 'centroid')
        about_code = 0;
    elseif strcmp(about, 'origin')
        about_code = 1;
    else
        error('mio:mass_properties', 'about must be ''centroid'' or ''origin''');
    end
    r = meshioplusplus_mex('mio_mass_properties', mesh, double(o.density), ...
        mio.internal.str(o.density_array, 'density_array'), ...
        mio.internal.str(o.region, 'region'), double(o.thickness), about_code);
    if isempty(r.warnings)
        r.warnings = {};
    else
        r.warnings = strsplit(r.warnings, char(10));
    end
end
