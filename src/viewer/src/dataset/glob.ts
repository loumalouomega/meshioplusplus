/**
 * A twin of the sequence glob language, for matching workspace file names on
 * the main thread — before any bytes exist to stage into the worker's MEMFS.
 *
 * The language is deliberately the narrow one `_sequence.py`/the C++ core
 * accept (doc/sequences.md): `*` matches any run of characters (possibly
 * empty) within a name, `?` exactly one character, everything else —
 * including `[`/`]` — is literal, there is no `**`, and the **directory part
 * of a pattern is literal** (globbing one directory at a time). The wasm
 * `sequenceEntries` remains the single authority on plan *ordering and
 * times* once files are staged; this twin only decides which workspace
 * files to stage, so a mismatch could over- or under-stage but never
 * mis-order a plan.
 */

/** Match a basename against a `*`/`?`-only pattern (classic two-pointer). */
export function nameMatch(pattern: string, name: string): boolean {
    let p = 0;
    let n = 0;
    let starP = -1;
    let starN = -1;
    while (n < name.length) {
        if (p < pattern.length && (pattern[p] === name[n] || pattern[p] === '?')) {
            p += 1;
            n += 1;
        } else if (p < pattern.length && pattern[p] === '*') {
            starP = p;
            starN = n;
            p += 1;
        } else if (starP >= 0) {
            p = starP + 1;
            starN += 1;
            n = starN;
        } else {
            return false;
        }
    }
    while (p < pattern.length && pattern[p] === '*') p += 1;
    return p === pattern.length;
}

/** Match a relative path against a pattern whose directory part is literal. */
export function globMatch(pattern: string, relPath: string): boolean {
    const slash = pattern.lastIndexOf('/');
    const dir = slash < 0 ? '' : pattern.slice(0, slash);
    const namePattern = slash < 0 ? pattern : pattern.slice(slash + 1);
    const pathSlash = relPath.lastIndexOf('/');
    const pathDir = pathSlash < 0 ? '' : relPath.slice(0, pathSlash);
    const name = pathSlash < 0 ? relPath : relPath.slice(pathSlash + 1);
    return dir === pathDir && nameMatch(namePattern, name);
}

export function isPattern(text: string): boolean {
    return text.includes('*') || text.includes('?');
}

/**
 * Natural-numeric comparison for display ordering: digit runs compare
 * numerically (by stripped length, then lexicographically — never parsed to
 * a number, so arbitrarily long runs are exact), other runs by char code, a
 * digit run before a non-digit run, ties broken on the raw strings. Display
 * only — the staged plan's order comes from the wasm side.
 */
export function naturalCompare(a: string, b: string): number {
    const runs = (s: string): string[] => s.match(/\d+|\D+/g) ?? [];
    const ra = runs(a);
    const rb = runs(b);
    for (let i = 0; i < Math.min(ra.length, rb.length); i++) {
        const x = ra[i];
        const y = rb[i];
        const xd = /^\d/.test(x);
        const yd = /^\d/.test(y);
        if (xd !== yd) return xd ? -1 : 1;
        if (xd) {
            const sx = x.replace(/^0+(?=\d)/, '');
            const sy = y.replace(/^0+(?=\d)/, '');
            if (sx.length !== sy.length) return sx.length - sy.length;
            if (sx !== sy) return sx < sy ? -1 : 1;
        } else if (x !== y) {
            return x < y ? -1 : 1;
        }
    }
    if (ra.length !== rb.length) return ra.length - rb.length;
    return a < b ? -1 : a > b ? 1 : 0;
}
