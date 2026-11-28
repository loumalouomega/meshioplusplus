function tf = has_boolean_provider()
%HAS_BOOLEAN_PROVIDER  Whether a surface boolean provider (Manifold) is registered.

    tf = meshioplusplus_mex('mio_has_boolean_provider');
end
