function tf = has_netgen()
%HAS_NETGEN  Whether the library carries the Netgen backend.
%   mio.tetrahedralize's 'netgen' backend needs a library built with
%   -DMESHIOPLUSPLUS_WITH_NETGEN=ON; this reports whether the loaded one was.

    tf = meshioplusplus_mex('mio_has_netgen');
end
