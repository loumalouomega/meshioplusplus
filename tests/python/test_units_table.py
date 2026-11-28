"""The Python twin's unit table (``_units._UNIT_TABLE``) against the core's
``field_unit_from_name`` / ``field_unit_from_text`` / ``field_unit_name``
(roadmap 9.10.1): the format readers and writers of both engines must turn the
same unit strings into the same ``meshio:unit:<array>`` entries."""

import numpy as np
import pytest

from meshioplusplus import _units

_core = pytest.importorskip("meshioplusplus._core")


@pytest.mark.parametrize("row", _units._UNIT_TABLE, ids=lambda r: r[0])
def test_every_python_row_matches_the_core(row):
    name, scale, exps = row
    np.testing.assert_array_equal(
        np.asarray(_core.field_unit_from_name(name)), [scale, *exps]
    )


@pytest.mark.parametrize(
    "text",
    ["m", " M/S ", "KG/M3", "kg/m3\x00\x00", "MPA", "mm/s^2", "furlong", "", "M2/S"],
)
def test_lenient_text_parsing_agrees(text):
    py = _units.field_unit_from_text(text)
    cpp = _core.field_unit_from_text(text)
    if cpp is None:
        assert py is None
    else:
        np.testing.assert_array_equal(py, np.asarray(cpp))


@pytest.mark.parametrize("row", _units._UNIT_TABLE, ids=lambda r: r[0])
def test_reverse_lookup_agrees(row):
    _, scale, exps = row
    unit = [scale, *exps]
    assert _units.field_unit_name(unit) == _core.field_unit_name(unit)


def test_reverse_lookup_of_an_unknown_unit_is_empty():
    unit = [7.0, 1, 0, 0, 0, 0, 0, 0]
    assert _units.field_unit_name(unit) == "" == _core.field_unit_name(unit)
