#!MC 1410
# Re-save a Tecplot SZL file as binary .plt, which meshio++ reads natively.
# Edit the two paths, then in Tecplot 360: Scripting > Play Macro/Script,
# or run `tec360 -b -p szplt_to_plt.mcr`.
# MIT License -- part of meshio++ (https://github.com/loumalouomega/meshioplusplus).
$!VARSET |IN| = 'run.szplt'
$!VARSET |OUT| = 'run.plt'
$!READDATASET  '"STANDARDSYNTAX" "1.0" "FILELIST_DATAFILES" "1" "|IN|"'
  DATASETREADER = 'Tecplot Subzone Data Loader'
  READDATAOPTION = NEW
  RESETSTYLE = YES
$!WRITEDATASET  "|OUT|"
  INCLUDETEXT = NO
  INCLUDEGEOM = NO
  INCLUDEDATASHARELINKAGE = YES
  BINARY = YES
  USEPOINTFORMAT = NO
  PRECISION = 9
  TECPLOTVERSIONTOWRITE = TECPLOTCURRENT
