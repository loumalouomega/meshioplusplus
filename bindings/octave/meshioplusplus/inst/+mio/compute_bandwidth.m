function bw = compute_bandwidth(mesh)
%COMPUTE_BANDWIDTH  The node-numbering bandwidth of a mesh (a number).

    bw = meshioplusplus_mex('mio_compute_bandwidth', mesh);
end
