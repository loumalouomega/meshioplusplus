"""The bundled single-file browser viewer.

``viewer.html`` is generated output, not source: it is the result of
``cd src/viewer && npm run build:embed``, copied here. See ``src/viewer/README.md``
for how to regenerate it and why it carries no WebAssembly.

This is a package rather than a bare data directory so that the file travels
with ``[tool.scikit-build] wheel.packages`` and linters do not treat the
directory as stray.
"""
