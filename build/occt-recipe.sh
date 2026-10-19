# occt-recipe.sh -- the one pinned OpenCASCADE source and configure recipe,
# sourced (not run) by build/build-occt.sh (the native static prefix the PyPI
# wheels and the release CLI binaries link) and by build/build-wasm-deps.sh
# --full (the wasm32 one), so the kernel every artifact ships is the same
# commit built from the same toolkit list.
#
# OCCT 7.6 is the floor the backend sets (Ubuntu 24.04's packages, the CI
# `occt` leg). It publishes no source tarball, so the source is always a
# shallow clone of the tag, refused unless it checks out exactly OCCT_COMMIT.

OCCT_VERSION="7.6.3"
OCCT_TAG="V7_6_3"
OCCT_COMMIT="b079fb9877ef64d4a8158a60fa157f59b096debb"
OCCT_URL="https://github.com/Open-Cascade-SAS/OCCT.git"

# Only the toolkits the core links (CMakeLists.txt's _mio_occt_libs): every
# module is off and these are named, which pulls in their dependency closure
# (modelling data and algorithms, shape healing, booleans, primitives) and
# nothing else -- no visualisation, OCAF, Draw, or third-party library.
OCCT_TOOLKITS="TKSTEP TKIGES TKMesh TKBool TKPrim"

# occt_source <download-dir> <source-dir>: the commit-verified clone of the
# tag (cached in <download-dir>) exported into <source-dir>, patched with
# build/patches/OCCT-<version>-*.patch.
occt_source() {
    _occt_git="$1/occt-$OCCT_VERSION.git"
    if [ ! -d "$_occt_git/.git" ] || [ "$(git -C "$_occt_git" rev-parse HEAD)" != "$OCCT_COMMIT" ]; then
        echo "  cloning: $OCCT_URL $OCCT_TAG"
        rm -rf "$_occt_git"
        git -c advice.detachedHead=false clone --quiet --depth 1 --branch "$OCCT_TAG" \
            "$OCCT_URL" "$_occt_git"
        _occt_got=$(git -C "$_occt_git" rev-parse HEAD)
        if [ "$_occt_got" != "$OCCT_COMMIT" ]; then
            rm -rf "$_occt_git"
            echo "error: $OCCT_URL $OCCT_TAG is $_occt_got, expected $OCCT_COMMIT" >&2
            exit 1
        fi
    else
        echo "  cached: occt-$OCCT_VERSION ($OCCT_TAG)"
    fi
    rm -rf "$2"
    mkdir -p "$2"
    (cd "$_occt_git" && git archive --format=tar HEAD) | tar -xf - -C "$2"
    for _occt_patch in "$OCCT_RECIPE_DIR/patches/OCCT-$OCCT_VERSION-"*.patch; do
        [ -e "$_occt_patch" ] || continue
        echo "  applying $(basename "$_occt_patch")"
        (cd "$2" && patch -p1 --forward --silent < "$_occt_patch")
    done
}

# occt_configure <source-dir> <build-dir> <prefix> <Static|Shared> [cmake args]:
# configure the toolkits above into <prefix> (Unix layout). $OCCT_CMAKE is the
# command run (default `cmake`; the wasm build sets `emcmake cmake`); the
# extra arguments come last, so they win.
#
# BUILD_RELEASE_DISABLE_EXCEPTIONS=OFF keeps OCCT's range and domain checks in
# a release build (upstream's default compiles them out with -DNo_Exception, as
# Debian's packages do not): fit_cad relies on a degenerate surface or curve
# throwing Standard_Failure, which demotes that one face, where the unchecked
# build reads out of bounds and crashes (BSplSLib::Resolution on a fillet's
# B-spline face, seen with 7.6.3).
occt_configure() {
    _occt_src="$1"; _occt_build="$2"; _occt_prefix="$3"; _occt_type="$4"
    shift 4
    ${OCCT_CMAKE:-cmake} -S "$_occt_src" -B "$_occt_build" \
        -DCMAKE_BUILD_TYPE=Release \
        -DINSTALL_DIR="$_occt_prefix" \
        -DINSTALL_DIR_LAYOUT=Unix \
        -DBUILD_LIBRARY_TYPE="$_occt_type" \
        -DBUILD_MODULE_FoundationClasses=OFF \
        -DBUILD_MODULE_ModelingData=OFF \
        -DBUILD_MODULE_ModelingAlgorithms=OFF \
        -DBUILD_MODULE_Visualization=OFF \
        -DBUILD_MODULE_ApplicationFramework=OFF \
        -DBUILD_MODULE_DataExchange=OFF \
        -DBUILD_MODULE_Draw=OFF \
        -DBUILD_ADDITIONAL_TOOLKITS="$OCCT_TOOLKITS" \
        -DBUILD_DOC_Overview=OFF \
        -DBUILD_RESOURCES=OFF \
        -DBUILD_USE_PCH=OFF \
        -DBUILD_RELEASE_DISABLE_EXCEPTIONS=OFF \
        -DUSE_FREETYPE=OFF -DUSE_FREEIMAGE=OFF -DUSE_OPENGL=OFF -DUSE_GLES2=OFF \
        -DUSE_TK=OFF -DUSE_TCL=OFF -DUSE_TBB=OFF -DUSE_VTK=OFF -DUSE_RAPIDJSON=OFF \
        -DUSE_DRACO=OFF -DUSE_FFMPEG=OFF -DUSE_OPENVR=OFF -DUSE_EIGEN=OFF \
        "$@"
}

# occt_licenses <source-dir> <prefix>: OCCT's licence texts into
# <prefix>/share/licenses/opencascade, for the notices a binary that links
# the kernel must carry (THIRD_PARTY_NOTICES.md).
occt_licenses() {
    mkdir -p "$2/share/licenses/opencascade"
    cp "$1/LICENSE_LGPL_21.txt" "$1/OCCT_LGPL_EXCEPTION.txt" "$2/share/licenses/opencascade/"
}
