function r = data_integrate(mesh, varargin)
%DATA_INTEGRATE  Cell-measure-weighted total and mean of cell data arrays.
%   R = mio.data_integrate(MESH) integrates every cell data array, weighting
%   each cell by its own length/area/volume -- mio.gradient's integration
%   counterpart. A cell whose measure is not computable (ragged, unsupported
%   type or degenerate), or whose value is non-finite in a component, is
%   excluded from that component's numerator and denominator. Reported for
%   the whole mesh (domain) and independently for every named Cell region
%   (regions): a cell in two regions contributes fully to both, one in none
%   to neither. A point-data-only name fails, naming point_data_to_cell_data
%   (mio.data_point_to_cell) as the fix. See doc/field_integration.md.
%   Option: 'names' (cell data arrays; [] for every one).
%
%   R is an n x 1 struct array with fields name, num_components, domain (a
%   struct: num_cells, num_skipped, components) and regions (an r x 1 struct
%   array: name, num_cells, num_skipped, components). Every components is a
%   (num_components x 4) matrix with columns total, mean, domain_measure and
%   num_nan.

    o = mio.internal.opts(varargin, {'names', []}, 'mio.data_integrate');
    r = meshioplusplus_mex('mio_data_integrate', mesh, mio.internal.cellstr_arg(o.names, 'names'));
end
