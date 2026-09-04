/**
 * Content hashing for manifest identity. The File System Access API never
 * reveals an absolute path, so a browser-side card can only be matched to a
 * server-side manifest (the companion process, doc/dashboard.md) by the
 * bytes both sides see. `crypto.subtle` exists only in secure contexts;
 * elsewhere the hash is simply unknown.
 */

export async function sha256Hex(text: string): Promise<string | null> {
    const subtle = globalThis.crypto?.subtle;
    if (!subtle) return null;
    try {
        const digest = await subtle.digest('SHA-256', new TextEncoder().encode(text));
        return [...new Uint8Array(digest)].map((b) => b.toString(16).padStart(2, '0')).join('');
    } catch {
        return null;
    }
}
