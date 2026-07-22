# Toolbar icons: tikz-ui/*.tex -> svg-ui/*.svg -> ../../viewer/src/ui/icons.ts
#
# Same visual language as the sibling CAD-Preview / MDPA-Preview icon sets, and
# the same toolchain the repository already uses for logo/ (pdflatex +
# dvisvgm). Note dvisvgm, not pdftocairo: see build-icons.mjs for why the
# post-processing differs from those projects'.
#
# `make ui` needs a TeX install; `make ts` does not, as long as svg-ui/*.svg is
# current -- which is why those are committed.
UI_TEXDIR   := tikz-ui
UI_BUILDDIR := build-ui
UI_SVGDIR   := svg-ui

UI_TEX := $(wildcard $(UI_TEXDIR)/*.tex)
UI_SVG := $(patsubst $(UI_TEXDIR)/%.tex,$(UI_SVGDIR)/%.svg,$(UI_TEX))

.PHONY: all ui ts clean
all: ts
ui: $(UI_SVG)

$(UI_SVGDIR)/%.svg: $(UI_TEXDIR)/%.tex preamble.tex
	@mkdir -p $(UI_BUILDDIR) $(UI_SVGDIR)
	cd $(UI_TEXDIR) && pdflatex -interaction=nonstopmode -halt-on-error \
	    -output-directory=../$(UI_BUILDDIR) $*.tex >/dev/null
	dvisvgm --pdf --no-fonts --output=$(UI_SVGDIR)/$*.svg $(UI_BUILDDIR)/$*.pdf >/dev/null 2>&1

# Regenerates ../viewer/src/ui/icons.ts. Pure Node, no LaTeX needed.
ts: ui
	node build-icons.mjs

clean:
	rm -rf $(UI_BUILDDIR)
