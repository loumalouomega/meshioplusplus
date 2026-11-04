function tf = has_mmg()
%HAS_MMG  Whether the library carries the MMG backend.
%   The MMG operations (mio.adapt, mio.discretize_levelset, mio.mesh_domain,
%   mio.move) need a library built with -DMESHIOPLUSPLUS_WITH_MMG=ON; this
%   reports whether the loaded one was. See also mio.has_mmg_elas.

    tf = meshioplusplus_mex('mio_has_mmg');
end
