function tf = has_crs_provider()
%HAS_CRS_PROVIDER  Whether a CRS provider (PROJ) is registered.

    tf = meshioplusplus_mex('mio_has_crs_provider');
end
