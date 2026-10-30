# The OpenCASCADE toolkits the OpenCASCADE provider links, as the core's own
# probe picks them (CMakeLists.txt). Sets _mio_occt_libs; OpenCASCADE must be found.
# The toolkits, as the core's own probe picks them: 7.8 renamed the STEP and
# IGES ones into TKDESTEP and TKDEIGES.
set(_mio_occt_libs TKernel TKMath TKG2d TKG3d TKBRep TKGeomBase TKGeomAlgo TKTopAlgo
                   TKShHealing TKBO TKBool TKMesh TKPrim)
if(TARGET TKDESTEP OR OpenCASCADE_VERSION VERSION_GREATER_EQUAL 7.8)
  list(APPEND _mio_occt_libs TKDESTEP TKXSBase TKDEIGES)
else()
  list(APPEND _mio_occt_libs TKSTEP TKSTEPBase TKSTEPAttr TKSTEP209 TKXSBase TKIGES)
endif()

