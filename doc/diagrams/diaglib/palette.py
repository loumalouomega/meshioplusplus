"""The one place colours, fonts and stroke widths live.

Values are the ``dataviz`` skill's brand-neutral reference palette (light
column). Categorical slots are assigned to *entities* and never cycled, so the
C++ core is the same blue in every figure and a refinement's green closure
cell is the same green on every page.
"""

# Surface and ink -----------------------------------------------------------
PAPER = "#fcfcfb"
PAPER_BORDER = "#dededa"
INK = "#0b0b0b"
INK_2 = "#52514e"
MUTED = "#898781"
HAIRLINE = "#e1e0d9"
WHITE = "#ffffff"

# Categorical slots, by entity ----------------------------------------------
CORE = "#2a78d6"  # slot 1 blue: the C++ core, C++ things, mid-edge nodes
PYTHON = "#eb6834"  # slot 2 orange: Python, face-centre nodes
FORMATS = "#1baf7a"  # slot 3 aqua: formats / the registry / files on disk
DATA = "#eda100"  # slot 4 yellow: data arrays
REGIONS = "#e87ba4"  # slot 5 magenta: regions / sets
WASM = "#008300"  # slot 6 green: WASM / browser / green closure cells
CABI = "#4a3aa7"  # slot 7 violet: the C ABI family (C, Fortran, Julia, R), body nodes
ERROR = "#e34948"  # slot 8 red: refusals, errors, red (fully split) cells

# Ordinal effort ramp (sequential blue) ------------------------------------
EFFORT = {"S": "#86b6ef", "M": "#3987e5", "L": "#1c5cab", "XL": "#0d366b"}

#: Every colour a figure may use; the tests refuse anything else.
ALL = frozenset(
    [
        PAPER,
        PAPER_BORDER,
        INK,
        INK_2,
        MUTED,
        HAIRLINE,
        WHITE,
        CORE,
        PYTHON,
        FORMATS,
        DATA,
        REGIONS,
        WASM,
        CABI,
        ERROR,
        *EFFORT.values(),
    ]
)

# Typography and strokes ---------------------------------------------------
FONT = "system-ui, -apple-system, 'Segoe UI', Roboto, sans-serif"
MONO = "ui-monospace, SFMono-Regular, Menlo, Consolas, 'Liberation Mono', monospace"
SIZE_LABEL = 13
SIZE_SMALL = 11
SIZE_TITLE = 15
STROKE_BOX = 1.5
STROKE_ARROW = 1.8
STROKE_HAIR = 1.0
FILL_OPACITY = 0.12
