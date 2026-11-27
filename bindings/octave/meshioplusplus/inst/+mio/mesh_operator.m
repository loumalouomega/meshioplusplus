function r = mesh_operator(mesh, varargin)
%MESH_OPERATOR  Assemble a P1 finite-element sparse CSR operator.
%   Options: 'kind' ('laplacian', 'mass_lumped', 'mass_consistent', 'gradient').
%   Returns a struct with rows, cols, nnz, indptr, indices, values.

    o = mio.internal.opts(varargin, {'kind', 'laplacian'}, 'mio.mesh_operator');
    switch o.kind
        case 'laplacian'
            k = 0;
        case 'mass_lumped'
            k = 1;
        case 'mass_consistent'
            k = 2;
        case 'gradient'
            k = 3;
        otherwise
            error('mio:mesh_operator', 'unknown kind ''%s''', o.kind);
    end
    r = meshioplusplus_mex('mio_mesh_operator', mesh, double(k));
end
