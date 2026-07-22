#!/bin/sh
# build.sh - regenerate the meshio++ logo assets from the TikZ sources.
#
#   ./build.sh
#
# Produces (committed): logo-with-text.svg, logo-icon.svg, logo.pdf,
# logo.png / logo-icon.png when a rasteriser is available, plus (via
# make_icon_assets.py, needs PyMuPDF+Pillow) logo-icon-square.png,
# social-preview.png, and doc/public/favicon.ico.
#
# Pipeline: gen_logo_tikz.py -> _mesh_icon.tex ; pdflatex -> PDF ;
# dvisvgm --pdf -> SVG ; PNG via PyMuPDF (fitz) / pdftoppm / convert if present.
set -eu

HERE=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
cd "$HERE"

PY="${PYTHON:-}"
if [ -z "$PY" ]; then
    if [ -x "$HERE/../../.venv/bin/python" ]; then
        PY="$HERE/../../.venv/bin/python"
    else
        PY=$(command -v python3 || command -v python)
    fi
fi

echo "== generating mesh-icon geometry =="
"$PY" gen_logo_tikz.py

compile_one() {
    tex="$1"      # logo | logo-icon
    echo "== pdflatex $tex.tex =="
    pdflatex -interaction=nonstopmode -halt-on-error "$tex.tex" >"$tex.build.log" 2>&1 \
        || { echo "pdflatex FAILED for $tex:"; tail -20 "$tex.build.log"; exit 1; }
    echo "== dvisvgm $tex.pdf -> svg =="
    dvisvgm --pdf --no-fonts --output="$tex.svg" "$tex.pdf" >>"$tex.build.log" 2>&1 \
        || { echo "dvisvgm FAILED for $tex:"; tail -20 "$tex.build.log"; exit 1; }
}

compile_one logo
compile_one logo-icon

# Canonical asset names.
cp -f logo.svg logo-with-text.svg

# PNG (best-effort). Prefer PyMuPDF (self-contained, no system libs).
rasterise() {
    src_pdf="$1"; dst_png="$2"
    if "$PY" - "$src_pdf" "$dst_png" <<'PYEOF' 2>/dev/null
import sys
try:
    import fitz  # PyMuPDF
except Exception:
    sys.exit(3)
doc = fitz.open(sys.argv[1])
pix = doc[0].get_pixmap(matrix=fitz.Matrix(4, 4), alpha=True)
pix.save(sys.argv[2])
PYEOF
    then
        echo "== PNG (PyMuPDF): $dst_png =="
    elif command -v pdftoppm >/dev/null 2>&1; then
        pdftoppm -png -r 300 -singlefile "$src_pdf" "${dst_png%.png}"
        echo "== PNG (pdftoppm): $dst_png =="
    elif command -v convert >/dev/null 2>&1; then
        convert -density 300 -background none "$src_pdf" "$dst_png"
        echo "== PNG (convert): $dst_png =="
    else
        echo "!! no PNG rasteriser (PyMuPDF/pdftoppm/convert) - skipping $dst_png"
    fi
}

rasterise logo.pdf logo.png
rasterise logo-icon.pdf logo-icon.png

# Derived small assets (square icon, favicon.ico, GitHub social-preview PNG).
# Best-effort: skip quietly if PyMuPDF/Pillow aren't available.
if "$PY" -c "import fitz, PIL" 2>/dev/null; then
    echo "== deriving icon assets (favicon.ico, logo-icon-square.png, social-preview.png) =="
    "$PY" make_icon_assets.py
else
    echo "!! PyMuPDF/Pillow not available - skipping favicon.ico/logo-icon-square.png/social-preview.png"
fi

# Tidy LaTeX aux (keep the committed .pdf/.svg/.png).
rm -f ./*.aux ./*.log ./*.build.log

echo
echo "== assets =="
ls -la logo-with-text.svg logo-icon.svg logo.pdf logo-icon.pdf 2>/dev/null || true
ls -la logo.png logo-icon.png logo-icon-square.png social-preview.png 2>/dev/null || true
ls -la ../public/favicon.ico 2>/dev/null || true
