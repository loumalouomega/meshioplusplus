"""``doc/roadmap.md`` checks itself: every item number cited in *Suggested
sequencing*, the dependency overview and the milestones table names an item
(a heading or a bold label) that exists, and every ``§N`` link names a section.

Only that planning text is checked. The history in *How this file works* and the
*Non-goals* cite items that closed on purpose, the open items quote files they
will create, and version numbers (``CMake 3.21``) look like item numbers.
"""

import pathlib
import re

REPO = pathlib.Path(__file__).resolve().parents[2]
ROADMAP = REPO / "doc" / "roadmap.md"
LINES = ROADMAP.read_text(encoding="utf-8").splitlines()

_HEADING = re.compile(r"^#{2,6}\s+(\d+(?:\.\d+)*)\.?\s")
_BOLD_ITEM = re.compile(r"\*\*(\d+(?:\.\d+)+)\b")
_CITATION = re.compile(r"(?<![\w.\-/§v])((?:[1-9]|1[0-4])\.\d+(?:\.\d+)*)(?![\w]|\.\d)")
_PATH = re.compile(
    r"`((?:src|tests|tools|bindings|doc|build|cmake|\.github|example|packages)/[\w./\-*<>{}]+)`"
)


def _item_ids():
    ids = set()
    for line in LINES:
        match = _HEADING.match(line)
        if match:
            ids.add(match.group(1))
        ids.update(_BOLD_ITEM.findall(line))
    return ids


def _section_span(start_prefix, end_prefix):
    start = next(i for i, line in enumerate(LINES) if line.startswith(start_prefix))
    end = next(
        i for i, line in enumerate(LINES) if i > start and line.startswith(end_prefix)
    )
    return start, end


def _planning_lines():
    """(1-based line number, text) of the text whose numbers must resolve."""
    sequencing, spec_index = _section_span(
        "## Suggested sequencing", "### Specification index"
    )
    for i in range(sequencing, spec_index):
        yield i + 1, LINES[i]


def _strip(line):
    line = re.sub(r"`[^`]*`", "", line)
    line = re.sub(r"\]\([^)]*\)", "]", line)
    return re.sub(r"https?://\S+", "", line)


def test_sequencing_and_dependencies_cite_existing_items():
    ids = _item_ids()
    stale = []
    for number, line in _planning_lines():
        if line.startswith("#"):
            continue
        for cited in _CITATION.findall(_strip(line)):
            parent = ".".join(cited.split(".")[:2])
            if cited not in ids and parent not in ids:
                stale.append((number, cited))
    assert stale == [], f"roadmap cites items that do not exist: {stale}"


def test_section_links_resolve():
    sections = {
        m.group(1) for line in LINES if (m := re.match(r"##\s+(\d+)\.\s", line))
    }
    cited = {m for line in LINES for m in re.findall(r"§(\d+)", line)}
    assert cited <= sections, cited - sections
