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


def _complete(script, words, cword, tmp_path):
    """COMPREPLY for `words` in a real bash with `script` sourced."""
    path = tmp_path / "c.bash"
    path.write_text(script)
    array = " ".join(f"'{w}'" for w in words)
    cmd = f"source {path}; COMP_WORDS=({array}); COMP_CWORD={cword}; _meshioplusplus; "
    cmd += 'printf "%s\\n" "${COMPREPLY[@]}"'
    r = subprocess.run(["bash", "-c", cmd], capture_output=True, text=True, check=True)
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


@pytest.mark.skipif(shutil.which("bash") is None, reason="needs bash")
@pytest.mark.parametrize("engine", ["python", "native"])
def test_bash_script_completes(engine, tmp_path):
    if engine == "native" and NATIVE is None:
        pytest.skip("no native CLI binary built")
    script = _python_script("bash") if engine == "python" else _native_script("bash")
    subprocess.run(["bash", "-n"], input=script, text=True, check=True)
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
        subprocess.run([shell, "-n"], input=script, text=True, check=True)


def test_the_native_verbs_options_match_their_parsers():
    """The native table comes from each verb's cli_parse call (describe mode)."""
    if NATIVE is None:
        pytest.skip("no native CLI binary built")
    script = _native_script("bash")
    line = next(ln for ln in script.splitlines() if ln.strip().startswith("geodesic)"))
    for opt in ("--source", "--source-region", "--time-factor", "--json", "-i", "-o"):
        assert f" {opt} " in line or f'"{opt} ' in line or f' {opt}"' in line
