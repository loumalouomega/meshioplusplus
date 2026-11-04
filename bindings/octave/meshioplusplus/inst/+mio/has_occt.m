function tf = has_occt()
%HAS_OCCT  Whether the loaded library carries the OpenCASCADE backend.
%   mio.fit_cad(..., 'backend', 'occt') needs a library built with
%   -DMESHIOPLUSPLUS_WITH_OCCT=ON. Returns true or false.

    tf = meshioplusplus_mex('mio_has_occt');
end
