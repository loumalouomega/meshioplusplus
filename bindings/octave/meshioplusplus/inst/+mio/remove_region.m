function mesh = remove_region(mesh, index)
%REMOVE_REGION  Remove the region at 1-based INDEX (in mio.regions order).

    meshioplusplus_mex('mio_remove_region', mesh, double(index));
end
