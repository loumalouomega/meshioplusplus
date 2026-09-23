"""The Fortran format-line parser shared by the fixed-column readers (the twin of
tests/cpp/test_keyword_card.cpp)."""

import pytest

from meshioplusplus import ReadError
from meshioplusplus.lsdyna._cards import parse_fortran_format, split_fixed


def _describe(fields):
    return " ".join(f"{k}{w}" for k, w in fields)


@pytest.mark.parametrize(
    "text, expected",
    [
        ("(3i8,6e20.13)", "i8 i8 i8 r20 r20 r20 r20 r20 r20"),
        ("(3i9,6e21.13e3)", "i9 i9 i9 r21 r21 r21 r21 r21 r21"),
        ("(1i7,2i9,6e21.13)", "i7 i9 i9 r21 r21 r21 r21 r21 r21"),
        ("(2i8,6g16.9)", "i8 i8 r16 r16 r16 r16 r16 r16"),
        ("(2(i8,e16.9))", "i8 r16 i8 r16"),
        ("(1P,3E20.12)", "r20 r20 r20"),
        ("( I5 , 2X , ES12.4 )", "i5 x2 r12"),
        ("i9", "i9"),
    ],
)
def test_expands_format_lines(text, expected):
    assert _describe(parse_fortran_format(text)) == expected


def test_etblock_and_repeats():
    fields = parse_fortran_format("(2i9,19a9)")
    assert len(fields) == 21 and fields[2] == ("a", 9)
    assert len(parse_fortran_format("(19i10)")) == 19


@pytest.mark.parametrize("bad", ["", "()", "(3q8)", "(i)", "(2i8", "(i8))"])
def test_refuses_garbage(bad):
    with pytest.raises(ReadError):
        parse_fortran_format(bad)


def test_splits_touching_fields():
    fields = parse_fortran_format("(3i8,6e20.13)")
    line = (
        "    5609       0       0 3.9797161316330E+00 2.5147820926190E-01"
        "-5.1500799817626E-01\r\n"
    )
    assert split_fixed(line, fields) == [
        "5609",
        "0",
        "0",
        "3.9797161316330E+00",
        "2.5147820926190E-01",
        "-5.1500799817626E-01",
    ]
    assert split_fixed("   12   3", parse_fortran_format("(i5,2x,i2)")) == ["12", "3"]
    assert split_fixed("", fields) == []
