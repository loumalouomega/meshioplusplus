function g = gmsh_info(info)
%GMSH_INFO  Copy a Gmsh format side channel.
%   G = mio.gmsh_info(INFO) returns a struct with
%     bounding_entities  n x 1 cell, the signed entity tags per cell block
%                        (Gmsh 4.1)
%     periodic           n x 1 struct array of links: entity (3 x 1:
%                        dimension, slave tag, master tag), affine (0 or 16
%                        coefficients, a column) and node_pairs (N x 2,
%                        1-based point rows: slave, master)
%   Entity tags remain file ids. Ordering and duplicates survive I/O;
%   operations do not remap these pairs. The copy remains valid after the
%   handle is released.

    g = meshioplusplus_mex('mio_gmsh_info', info);
end
