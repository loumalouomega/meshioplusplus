"""Tests for ``tools/check-abi-version.sh``.

This script is the *only* gate for Tier B ABI breaks -- an edit to the body of an
existing inline function -- because ``tests/cpp/test_abi_layout.cpp`` provably
cannot see one: a ``sizeof``/``alignof`` snapshot is unchanged by an inline-body
edit. So it is exactly the kind of check whose failure mode is silence.

It shipped that way once. Through v9.4.0 the review lookup matched any row of
``doc/abi_reviews.md`` whose first column was the current ABI number, and two
ABI-3 rows already existed -- so every later header change matched one of them
and exited 0 announcing "records the additive review", naming a review of a
different release entirely. Nothing tested that the gate *fired*, which is the
same lesson ``detail/mesh_backend_check.hpp`` learned in v9.1.0 when its guard
shipped inert for a full release.

Each test therefore drives the real script over a throwaway git repository and
asserts the exit status **and** the diagnostic. Nothing here touches the real
tree (the ``tests/python/test_gen_doc_images.py`` convention).
"""

import pathlib
import shutil
import subprocess

import pytest

REPO = pathlib.Path(__file__).resolve().parent.parent.parent
SCRIPT = REPO / "tools" / "check-abi-version.sh"


def _usable_bash():
    """Is there a ``bash`` on PATH that can actually run a POSIX script?

    Probed rather than inferred from ``sys.platform``. On GitHub's
    ``windows-latest``, ``bash`` resolves to ``C:\\Windows\\System32\\bash.exe``
    -- the WSL launcher -- ahead of Git Bash, and with no distribution installed
    it fails with a UTF-16 "Windows Subsystem for Linux has no installed
    distributions." A Windows box whose PATH does reach Git Bash runs these
    tests fine, which a platform check would have skipped for no reason.

    Nothing is lost where this skips: the gate runs in the Linux ``lint`` job.
    """
    try:
        probe = subprocess.run(
            ["bash", "-c", "printf ok"],
            capture_output=True,
            text=True,
            timeout=60,
        )
    except (OSError, subprocess.SubprocessError):
        return False
    return probe.returncode == 0 and probe.stdout.strip() == "ok"


pytestmark = pytest.mark.skipif(
    not _usable_bash(),
    reason="no usable POSIX bash on PATH (e.g. the Windows WSL stub)",
)

BASE_VERSION = "9.4.0"
NEXT_VERSION = "9.4.1"

# The fixture repo mirrors only what the script reads: the ABI header, the
# project VERSION line, the reviews table, and some installed headers.
CMAKELISTS = """cmake_minimum_required(VERSION 3.15...3.30)

project(
  meshioplusplus_core
  VERSION {version}
  LANGUAGES C CXX)
"""

ABI_HEADER = "#define MESHIOPLUSPLUS_ABI_VERSION {abi}\n"

REVIEWS = """# ABI additive-change reviews

| ABI | releases | headers changed | why it is additive |
| --- | --- | --- | --- |
{rows}"""

ROW = "| {abi} | v{version} | {headers} | Reviewed as additive. |\n"


def _git(repo, *args):
    subprocess.run(
        ["git", *args],
        cwd=repo,
        check=True,
        capture_output=True,
        text=True,
    )


def _write(repo, relpath, text):
    path = repo / relpath
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


def _set_version(repo, version):
    _write(repo, "CMakeLists.txt", CMAKELISTS.format(version=version))


def _set_abi(repo, abi):
    _write(
        repo,
        "src/cpp/include/meshioplusplus/abi_version.hpp",
        ABI_HEADER.format(abi=abi),
    )


def _set_reviews(repo, rows):
    _write(repo, "doc/abi_reviews.md", REVIEWS.format(rows="".join(rows)))


def _run(repo, base=None):
    cmd = ["bash", str(repo / "tools" / "check-abi-version.sh")]
    if base is not None:
        cmd.append(base)
    return subprocess.run(cmd, cwd=repo, capture_output=True, text=True)


@pytest.fixture
def repo(tmp_path):
    """A tagged v9.4.0 baseline: ABI 3, two installed headers, one review row.

    The review row deliberately covers the *previous* release only -- that is the
    shape that used to pass vacuously.
    """
    root = tmp_path / "fixture"
    root.mkdir()
    shutil.copytree(SCRIPT.parent, root / "tools", dirs_exist_ok=True)

    _git(root, "init", "-q", "-b", "main")
    _git(root, "config", "user.email", "test@example.com")
    _git(root, "config", "user.name", "Test")

    _set_version(root, BASE_VERSION)
    _set_abi(root, 3)
    _set_reviews(root, [ROW.format(abi=3, version=BASE_VERSION, headers="`mesh.hpp`")])
    _write(root, "src/cpp/include/meshioplusplus/mesh.hpp", "// mesh\n")
    _write(root, "src/cpp/include/meshioplusplus/region.hpp", "// region\n")

    _git(root, "add", "-A")
    _git(root, "commit", "-qm", "baseline")
    _git(root, "tag", f"v{BASE_VERSION}")
    return root


def test_unreviewed_header_change_fails(repo):
    """The reported bug: a header changes, the ABI holds, and the only review on
    file is for the *previous* release. This used to exit 0."""
    _write(repo, "src/cpp/include/meshioplusplus/region.hpp", "// region, edited\n")
    _set_version(repo, NEXT_VERSION)

    result = _run(repo, f"v{BASE_VERSION}")

    assert result.returncode == 1, result.stdout + result.stderr
    assert "records no" in result.stderr
    assert f"v{NEXT_VERSION}" in result.stderr
    # The point of the fix: an existing row for a different release must not count.
    assert "some *other* release does not count" in result.stderr


def test_review_row_for_previous_release_does_not_cover_a_new_change(repo):
    """The same shape without a version bump -- version_now stays at the previous
    release, whose row exists. Keying on the release alone would pass here; the
    header-coverage half is what fails it."""
    _write(repo, "src/cpp/include/meshioplusplus/region.hpp", "// region, edited\n")

    result = _run(repo, f"v{BASE_VERSION}")

    assert result.returncode == 1, result.stdout + result.stderr
    assert "does not" in result.stderr and "account for" in result.stderr
    assert "region.hpp" in result.stderr


def test_reviewed_header_change_passes(repo):
    """A row for this release naming this header is the supported way through."""
    _write(repo, "src/cpp/include/meshioplusplus/region.hpp", "// region, edited\n")
    _set_version(repo, NEXT_VERSION)
    _set_reviews(
        repo,
        [
            ROW.format(abi=3, version=BASE_VERSION, headers="`mesh.hpp`"),
            ROW.format(abi=3, version=NEXT_VERSION, headers="`region.hpp`"),
        ],
    )

    result = _run(repo, f"v{BASE_VERSION}")

    assert result.returncode == 0, result.stdout + result.stderr
    assert "reviews them for" in result.stdout


def test_partial_review_fails(repo):
    """A row that names one changed header but not the other is not a review of
    this change. Catches a second edit riding an earlier row."""
    _write(repo, "src/cpp/include/meshioplusplus/region.hpp", "// region, edited\n")
    _write(repo, "src/cpp/include/meshioplusplus/mesh.hpp", "// mesh, edited\n")
    _set_version(repo, NEXT_VERSION)
    _set_reviews(
        repo,
        [
            ROW.format(abi=3, version=BASE_VERSION, headers="`mesh.hpp`"),
            ROW.format(abi=3, version=NEXT_VERSION, headers="`region.hpp`"),
        ],
    )

    result = _run(repo, f"v{BASE_VERSION}")

    assert result.returncode == 1, result.stdout + result.stderr
    assert "mesh.hpp" in result.stderr
    # region.hpp IS reviewed, so it must not be listed as missing.
    missing = result.stderr.split("account for every installed header")[1]
    assert "region.hpp" not in missing.split("Either name them")[0]


def test_abi_bump_passes_without_any_review(repo):
    """Bumping the ABI is the other legitimate answer, and needs no review row."""
    _write(repo, "src/cpp/include/meshioplusplus/region.hpp", "// region, edited\n")
    _set_abi(repo, 4)
    _set_version(repo, NEXT_VERSION)

    result = _run(repo, f"v{BASE_VERSION}")

    assert result.returncode == 0, result.stdout + result.stderr
    assert "ABI version was bumped" in result.stdout


def test_no_header_change_passes(repo):
    """Changing anything outside the installed headers is not this gate's business."""
    _write(repo, "README.md", "docs\n")
    _set_version(repo, NEXT_VERSION)

    result = _run(repo, f"v{BASE_VERSION}")

    assert result.returncode == 0, result.stdout + result.stderr
    assert "no installed header changed" in result.stdout


def test_unparseable_version_is_an_error_not_a_pass(repo):
    """A gate that cannot read the version must fail loudly rather than skip."""
    _write(repo, "src/cpp/include/meshioplusplus/region.hpp", "// region, edited\n")
    _write(repo, "CMakeLists.txt", "project(meshioplusplus_core LANGUAGES CXX)\n")

    result = _run(repo, f"v{BASE_VERSION}")

    assert result.returncode == 1, result.stdout + result.stderr
    assert "cannot parse the project VERSION" in result.stderr
