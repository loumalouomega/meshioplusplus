"""The ``completion`` verb: a bash, zsh or fish completion script for this CLI,
generated from its own argparse parser, so the completions follow the verbs
and their options as they change (roadmap 7.6.2). The native CLI's
``completion`` verb writes the same scripts from its own parsers.

Install, for example::

    meshioplusplus completion bash > ~/.local/share/bash-completion/completions/meshioplusplus
    meshioplusplus completion zsh > "${fpath[1]}/_meshioplusplus"
    meshioplusplus completion fish > ~/.config/fish/completions/meshioplusplus.fish
"""

import argparse

PROG = "meshioplusplus"


def add_args(parser):
    parser.add_argument("shell", choices=("bash", "zsh", "fish"), help="the shell")


def _subparsers(parser):
    for action in parser._actions:
        if isinstance(action, argparse._SubParsersAction):
            return action
    return None


def _options(parser):
    out = set()
    for action in parser._actions:
        if action.option_strings and not isinstance(action, argparse._HelpAction):
            out.update(action.option_strings)
    return sorted(out)


def verbs(parser):
    """``[(name, aliases, options, subverbs)]`` for every verb, in order."""
    action = _subparsers(parser)
    seen = {}
    order = []
    for name, sub in action.choices.items():
        key = id(sub)
        if key in seen:
            seen[key][1].append(name)  # an alias of a verb already listed
            continue
        options, subverbs = _options(sub), []
        # A first positional with fixed choices (`completion bash|zsh|fish`)
        # completes like a sub-verb.
        positional = [a for a in sub._actions if not a.option_strings]
        if (
            positional
            and positional[0].choices
            and not isinstance(positional[0], argparse._SubParsersAction)
        ):
            subverbs = [str(c) for c in positional[0].choices]
        inner = _subparsers(sub)
        if inner is not None:
            for sub_name, sub_parser in inner.choices.items():
                subverbs.append(sub_name)
                options = sorted(set(options) | set(_options(sub_parser)))
        entry = (name, [], options, subverbs)
        seen[key] = entry
        order.append(entry)
    return order


def script(shell, entries, prog=PROG):
    """The completion script for ``shell`` over ``verbs(parser)``'s entries."""
    names = []
    for name, aliases, _, _ in entries:
        names += [name] + aliases
    lines = []
    if shell == "bash":
        lines += [
            f"# bash completion for {prog} ({prog} completion bash)",
            "_meshioplusplus() {",
            '    local cur="${COMP_WORDS[COMP_CWORD]}" opts="" subs=""',
            '    if [ "$COMP_CWORD" -eq 1 ]; then',
            f'        COMPREPLY=($(compgen -W "{" ".join(names)}" -- "$cur"))',
            "        return",
            "    fi",
            '    case "${COMP_WORDS[1]}" in',
        ]
        for name, aliases, options, subverbs in entries:
            pattern = "|".join([name] + aliases)
            lines.append(
                f'        {pattern}) opts="{" ".join(options)}" subs="{" ".join(subverbs)}";;'
            )
        lines += [
            "    esac",
            '    if [ "$COMP_CWORD" -eq 2 ] && [ -n "$subs" ]; then',
            '        COMPREPLY=($(compgen -W "$subs" -- "$cur"))',
            '    elif [[ "$cur" == -* ]]; then',
            '        COMPREPLY=($(compgen -W "$opts" -- "$cur"))',
            "    else",
            '        COMPREPLY=($(compgen -f -- "$cur"))',
            "    fi",
            "}",
            f"complete -o filenames -F _meshioplusplus {prog}",
        ]
    elif shell == "zsh":
        lines += [
            f"#compdef {prog}",
            f"# zsh completion for {prog} ({prog} completion zsh)",
            "_meshioplusplus() {",
            "    local -a opts subs",
            "    if (( CURRENT == 2 )); then",
            f"        compadd -- {' '.join(names)}",
            "        return",
            "    fi",
            "    case $words[2] in",
        ]
        for name, aliases, options, subverbs in entries:
            pattern = "|".join([name] + aliases)
            lines.append(
                f"        {pattern}) opts=({' '.join(options)}) subs=({' '.join(subverbs)});;"
            )
        lines += [
            "    esac",
            "    if (( CURRENT == 3 )) && (( ${#subs} )); then",
            "        compadd -- $subs",
            "    elif [[ $PREFIX == -* ]]; then",
            "        compadd -- $opts",
            "    else",
            "        _files",
            "    fi",
            "}",
            f"compdef _meshioplusplus {prog}",
        ]
    else:
        lines += [
            f"# fish completion for {prog} ({prog} completion fish)",
            f"complete -c {prog} -f",
            f"complete -c {prog} -n __fish_use_subcommand -a '{' '.join(names)}'",
        ]
        for name, aliases, options, subverbs in entries:
            seen = " ".join([name] + aliases)
            cond = f"-n '__fish_seen_subcommand_from {seen}'"
            if subverbs:
                lines.append(f"complete -c {prog} {cond} -a '{' '.join(subverbs)}'")
            for opt in options:
                if opt.startswith("--"):
                    flag = f"-l {opt[2:]}"
                elif len(opt) == 2:
                    flag = f"-s {opt[1:]}"
                else:
                    flag = f"-o {opt[1:]}"
                lines.append(f"complete -c {prog} {cond} {flag}")
            lines.append(f"complete -c {prog} {cond} -F")
    return "\n".join(lines) + "\n"


def completion_cmd(args, parser):
    print(script(args.shell, verbs(parser)), end="")
    return 0
