/**
 * Source-spec suggestion for the add-case flow: given the workspace-relative
 * paths a user selected, propose the tightest `SourceSpec`.
 *
 * - one file — `Path`;
 * - several files in one directory whose names share a prefix and suffix
 *   around a varying run — a `Pattern` (`prefix*suffix`), but only when the
 *   pattern matches the selection *exactly* against the full listing (a
 *   pattern silently matching future files is a feature; one matching
 *   *other existing* files is a wrong suggestion);
 * - anything else — an explicit `Paths` list in natural-numeric order.
 *
 * Pure string work, unit-tested; the user can always override the result.
 */

import { globMatch, naturalCompare } from './glob';
import type { SourceSpec } from './manifest';

function dirOf(path: string): string {
    const slash = path.lastIndexOf('/');
    return slash < 0 ? '' : path.slice(0, slash);
}

function nameOf(path: string): string {
    const slash = path.lastIndexOf('/');
    return slash < 0 ? path : path.slice(slash + 1);
}

function commonPrefix(names: string[]): string {
    let prefix = names[0];
    for (const name of names) {
        let i = 0;
        while (i < prefix.length && i < name.length && prefix[i] === name[i]) i += 1;
        prefix = prefix.slice(0, i);
    }
    return prefix;
}

function commonSuffix(names: string[]): string {
    let suffix = names[0];
    for (const name of names) {
        let i = 0;
        while (
            i < suffix.length &&
            i < name.length &&
            suffix[suffix.length - 1 - i] === name[name.length - 1 - i]
        ) {
            i += 1;
        }
        suffix = suffix.slice(suffix.length - i);
    }
    return suffix;
}

/**
 * Suggest a source for `selected`, verified against `allFiles` (the whole
 * workspace listing) so a suggested pattern never matches anything the user
 * did not pick.
 */
export function suggestSource(selected: string[], allFiles: string[]): SourceSpec {
    if (selected.length === 0) {
        throw new Error('meshio++: dataset: nothing selected');
    }
    const sorted = [...selected].sort(naturalCompare);
    if (sorted.length === 1) return { Path: sorted[0] };

    const dirs = new Set(sorted.map(dirOf));
    if (dirs.size === 1) {
        const names = sorted.map(nameOf);
        const prefix = commonPrefix(names);
        const suffix = commonSuffix(names);
        // The varying middle must be non-empty for at least one name and the
        // prefix/suffix must not overlap on the shortest name.
        const shortest = Math.min(...names.map((n) => n.length));
        if (prefix.length + suffix.length <= shortest) {
            const dir = dirs.values().next().value ?? '';
            const pattern = (dir ? `${dir}/` : '') + `${prefix}*${suffix}`;
            const matched = allFiles.filter((f) => globMatch(pattern, f));
            const wanted = new Set(sorted);
            if (matched.length === wanted.size && matched.every((f) => wanted.has(f))) {
                return { Pattern: pattern };
            }
        }
    }
    return { Paths: sorted };
}
