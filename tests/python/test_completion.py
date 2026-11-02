"""Shell completions for both CLIs (roadmap 7.6.2), generated from their own
parsers: the scripts parse, list every verb with its options, and complete in
a real bash."""

import contextlib
import io
import shutil
import subprocess

import pytest

import meshioplusplus as mp
from meshioplusplus._cli import _completion
from meshioplusplus._cli._main import build_parser

from .test_cli_json import NATIVE


def _python_script(shell):
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        assert mp._cli.main(["completion", shell]) in (0, None)
    return out.getvalue()


def _native_script(shell):
    r = subprocess.run([NATIVE, "completion", shell], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    return r.stdout


# The bash the tests drive, by its full path: on Windows a bare `bash` is looked
# up in System32 before PATH, and finds the WSL launcher there.
BASH = shutil.which("bash")


def _complete(script, words, cword, tmp_path):
    """COMPREPLY for `words` in a real bash with `script` sourced."""
    path = tmp_path / "c.bash"
    path.write_bytes(script.encode())  # no CRLF translation on Windows
    array = " ".join(f"'{w}'" for w in words)
    cmd = f"source '{path.as_posix()}'; COMP_WORDS=({array}); COMP_CWORD={cword}; _meshioplusplus; "
    cmd += 'printf "%s\\n" "${COMPREPLY[@]}"'
    r = subprocess.run([BASH, "-c", cmd], capture_output=True, text=True, check=True)
    return r.stdout.split()


def test_every_python_verb_and_option_is_listed():
    entries = _completion.verbs(build_parser())
    names = {e[0] for e in entries}
    assert {"convert", "graph", "geodesic", "data", "completion"} <= names
    graph = next(e for e in entries if e[0] == "graph")
    assert {"--kind", "--format", "--json"} <= set(graph[2])
    convert = next(e for e in entries if e[0] == "convert")
    assert "c" in convert[1]  # an alias
    data = next(e for e in entries if e[0] == "data")
    assert "info" in data[3]
    assert next(e for e in entries if e[0] == "completion")[3] == [
        "bash",
        "zsh",
        "fish",
    ]


def _bash_usable():
    """A bash that runs a command. On Windows the first `bash` on PATH can be
    the WSL launcher, which with no distribution installed prints a (UTF-16)
    message and still exits 0 for `-c true`, so the output is checked too."""
    if BASH is None:
        return False
    try:
        r = subprocess.run([BASH, "-c", "echo ok"], capture_output=True)
    except OSError:
        return False
    return r.returncode == 0 and r.stdout.strip() == b"ok"


@pytest.mark.skipif(not _bash_usable(), reason="needs a working bash")
@pytest.mark.parametrize("engine", ["python", "native"])
def test_bash_script_completes(engine, tmp_path):
    if engine == "native" and NATIVE is None:
        pytest.skip("no native CLI binary built")
    script = _python_script("bash") if engine == "python" else _native_script("bash")
    subprocess.run([BASH, "-n"], input=script.encode(), check=True)
    assert _complete(script, ["meshioplusplus", "geod"], 1, tmp_path) == ["geodesic"]
    assert _complete(script, ["meshioplusplus", "graph", "--ki"], 2, tmp_path) == [
        "--kind"
    ]
    assert "gradient" in _complete(
        script, ["meshioplusplus", "data", "gr"], 2, tmp_path
    )
    assert _complete(script, ["meshioplusplus", "completion", "z"], 2, tmp_path) == [
        "zsh"
    ]


@pytest.mark.parametrize("shell", ["zsh", "fish"])
@pytest.mark.parametrize("engine", ["python", "native"])
def test_zsh_and_fish_scripts(engine, shell):
    if engine == "native" and NATIVE is None:
        pytest.skip("no native CLI binary built")
    script = _python_script(shell) if engine == "python" else _native_script(shell)
    if shell == "zsh":
        assert script.startswith("#compdef meshioplusplus")
        assert "        graph) opts=(" in script
    else:
        assert (
            "complete -c meshioplusplus -n '__fish_seen_subcommand_from graph' -l kind"
            in script
        )
    if shutil.which(shell):
        subprocess.run([shutil.which(shell), "-n"], input=script.encode(), check=True)


def test_the_native_verbs_options_match_their_parsers():
    """The native table comes from each verb's cli_parse call (describe mode)."""
    if NATIVE is None:
        pytest.skip("no native CLI binary built")
    script = _native_script("bash")
    line = next(ln for ln in script.splitlines() if ln.strip().startswith("geodesic)"))
    for opt in ("--source", "--source-region", "--time-factor", "--json", "-i", "-o"):
        assert f" {opt} " in line or f'"{opt} ' in line or f' {opt}"' in line
