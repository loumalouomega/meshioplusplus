@echo off
rem configure.bat - configure (and optionally build) the meshio C++ core
rem standalone on Windows (MSVC). Run from anywhere; the build tree is created
rem next to this script (build\cpp-<build-type>).
rem
rem   configure.bat
rem   configure.bat --backend OPENMP --tests --build
rem
rem Notes:
rem   * The STL backend (default) needs nothing extra on MSVC - its parallel
rem     algorithms are built into the standard library.
rem   * HDF5/netCDF/zlib default to OFF on Windows (the Python fallbacks and
rem     the h5py/netCDF4 wheels cover those formats); pass --with-hdf5 etc.
rem     if you have the development libraries installed.

setlocal EnableDelayedExpansion

set "SCRIPT_DIR=%~dp0"
for %%I in ("%SCRIPT_DIR%..") do set "SOURCE_DIR=%%~fI"

set "BACKEND=STL"
set "BUILD_TYPE=Release"
set "WITH_HDF5=OFF"
set "WITH_NETCDF=OFF"
set "WITH_ZLIB=OFF"
set "TESTS=OFF"
set "DO_BUILD=no"
set "PYTHON_EXE="
set "TBB_DIR="

:parse
if "%~1"=="" goto endparse
if /I "%~1"=="--backend"        set "BACKEND=%~2" & shift & shift & goto parse
if /I "%~1"=="--build-type"     set "BUILD_TYPE=%~2" & shift & shift & goto parse
if /I "%~1"=="--with-hdf5"      set "WITH_HDF5=ON" & shift & goto parse
if /I "%~1"=="--without-hdf5"   set "WITH_HDF5=OFF" & shift & goto parse
if /I "%~1"=="--with-netcdf"    set "WITH_NETCDF=ON" & shift & goto parse
if /I "%~1"=="--without-netcdf" set "WITH_NETCDF=OFF" & shift & goto parse
if /I "%~1"=="--with-zlib"      set "WITH_ZLIB=ON" & shift & goto parse
if /I "%~1"=="--without-zlib"   set "WITH_ZLIB=OFF" & shift & goto parse
if /I "%~1"=="--tests"          set "TESTS=ON" & shift & goto parse
if /I "%~1"=="--build"          set "DO_BUILD=yes" & shift & goto parse
if /I "%~1"=="--python"         set "PYTHON_EXE=%~2" & shift & shift & goto parse
if /I "%~1"=="--tbb-dir"        set "TBB_DIR=%~2" & shift & shift & goto parse
if /I "%~1"=="-h" goto usage
if /I "%~1"=="--help" goto usage
echo Unknown option: %~1
goto usage

:endparse

if "%PYTHON_EXE%"=="" (
    if exist "%SOURCE_DIR%\.venv\Scripts\python.exe" (
        set "PYTHON_EXE=%SOURCE_DIR%\.venv\Scripts\python.exe"
    ) else (
        set "PYTHON_EXE=python"
    )
)

set "BUILD_DIR=%SCRIPT_DIR%cpp-%BUILD_TYPE%"

set "EXTRA="
for /f "delims=" %%P in ('"%PYTHON_EXE%" -c "import pybind11; print(pybind11.get_cmake_dir())" 2^>nul') do set "EXTRA=-Dpybind11_DIR=%%P"
if not "%TBB_DIR%"=="" set "EXTRA=%EXTRA% -DTBB_DIR=%TBB_DIR%"

echo == meshio configure ==
echo   source:    %SOURCE_DIR%
echo   build:     %BUILD_DIR%
echo   type:      %BUILD_TYPE%
echo   backend:   %BACKEND%
echo   HDF5:      %WITH_HDF5%   netCDF: %WITH_NETCDF%   zlib: %WITH_ZLIB%
echo   tests:     %TESTS%
echo   python:    %PYTHON_EXE%
echo.

cmake -S "%SOURCE_DIR%" -B "%BUILD_DIR%" ^
    -DMESHIO_PARALLEL_BACKEND=%BACKEND% ^
    -DMESHIO_WITH_HDF5=%WITH_HDF5% ^
    -DMESHIO_WITH_NETCDF=%WITH_NETCDF% ^
    -DMESHIO_WITH_ZLIB=%WITH_ZLIB% ^
    -DMESHIO_BUILD_TESTS=%TESTS% ^
    -DPython_EXECUTABLE="%PYTHON_EXE%" ^
    %EXTRA%
if errorlevel 1 exit /b 1

echo.
echo == next steps ==
echo   cmake --build "%BUILD_DIR%" --config %BUILD_TYPE%
if "%TESTS%"=="ON" echo   ctest --test-dir "%BUILD_DIR%" -C %BUILD_TYPE% --output-on-failure
echo.
echo Python package (editable) with the same options:
echo   set CMAKE_ARGS=-DMESHIO_PARALLEL_BACKEND=%BACKEND% -DMESHIO_WITH_HDF5=%WITH_HDF5% -DMESHIO_WITH_NETCDF=%WITH_NETCDF% -DMESHIO_WITH_ZLIB=%WITH_ZLIB%
echo   pip install --no-build-isolation -e "%SOURCE_DIR%"

if "%DO_BUILD%"=="yes" (
    echo.
    echo == building ==
    cmake --build "%BUILD_DIR%" --config %BUILD_TYPE%
)
exit /b 0

:usage
echo Usage: configure.bat [options]
echo   --backend ^<SEQ^|STL^|OPENMP^|TBB^>  parallel backend (default: STL)
echo   --build-type ^<type^>             CMake build type (default: Release)
echo   --with-hdf5 / --without-hdf5     HDF5-backed formats (default: off on Windows)
echo   --with-netcdf / --without-netcdf
echo   --with-zlib / --without-zlib
echo   --tests                          also build the GoogleTest suite (CTest)
echo   --build                          run the build after configuring
echo   --python ^<exe^>                  Python executable (default: auto)
echo   --tbb-dir ^<path^>                TBBConfig.cmake dir
exit /b 1
