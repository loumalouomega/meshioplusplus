"""Shared machinery for the generated documentation diagrams.

Everything under ``doc/public/diagrams/`` is emitted by ``gen_diagrams.py``
through this package: one palette, one SVG builder, one projector and one set
of code-table loaders, so every figure shares a look and can be regenerated
byte for byte. Standard library only -- no numpy, no meshioplusplus import --
so it runs wherever CI runs the ``lint`` job.
"""
