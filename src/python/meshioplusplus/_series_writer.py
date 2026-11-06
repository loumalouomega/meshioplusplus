"""The shared step/close/context-manager half of the core-backed ``SeriesWriter`` shims.

A format's ``SeriesWriter`` opens either the C++ core's series writer
(``self._core``) or, when the core declines (``_fallback.core_declined``), the
Python twin (``self._py``) in its own ``__init__``; this base forwards each step
to whichever is open and finalizes the core file on close when it holds a step.
"""


class CoreSeriesWriter:
    def write(self, time, mesh):
        if self._core is not None:
            self._core.write(float(time), mesh)
        else:
            self._py.write(time, mesh)

    def close(self):
        if self._core is not None:
            if self._core.num_steps():
                self._core.finalize()
            self._core = None
        if self._py is not None:
            self._py.close()
            self._py = None

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
        return False
