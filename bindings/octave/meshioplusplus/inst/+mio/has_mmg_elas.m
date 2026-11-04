function tf = has_mmg_elas()
%HAS_MMG_ELAS  Whether the library's MMG has ELAS, which mio.move needs.

    tf = meshioplusplus_mex('mio_has_mmg_elas');
end
