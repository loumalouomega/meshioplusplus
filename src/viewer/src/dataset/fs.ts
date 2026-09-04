/**
 * The workspace: how the dataset page reaches the user's local files.
 *
 * One interface, two implementations, chosen by capability at pick time:
 *
 * - **File System Access** (Chromium/Edge): `showDirectoryPicker()` gives a
 *   real directory handle — files enumerate lazily, and the manifest saves
 *   **back in place** through `createWritable()`, which is the roadmap's
 *   "mix hand edits and UI edits on the same file" story. Permission is
 *   re-requested on save when a reload dropped it, and the handle is
 *   persisted (best-effort, `persist.ts`) so the next visit offers a
 *   "Reopen" button — restore MUST run behind that click, because a
 *   persisted handle comes back with permission `'prompt'` and
 *   `requestPermission` demands a user gesture.
 * - **Fallback** (`<input type="file" webkitdirectory>`): every browser can
 *   *read* a directory this way (`webkitRelativePath` gives the same
 *   relative paths), but nothing can write back — saving downloads a copy,
 *   stated in the UI and documented in doc/datasets.md.
 *
 * The picked directory IS the manifest's directory (`base_dir`): entry
 * sources are workspace-relative, so Python's `DatasetManifest.load`
 * resolves them identically (doc/datasets.md's relative-path rule).
 */

export interface WorkspaceFile {
    relPath: string;
    bytes(): Promise<ArrayBuffer>;
    /** Epoch milliseconds of the file's last modification, or null when the
     * backend cannot say (the overview cards and the scan cache key on it). */
    lastModified(): Promise<number | null>;
}

export interface Workspace {
    kind: 'fsa' | 'fallback';
    /** The picked directory handle (FSA only) — what `persist.ts` stores. */
    readonly root?: FileSystemDirectoryHandle;
    /** Every file under the picked directory, workspace-relative. */
    files: WorkspaceFile[];
    /** `*.json` candidates at the workspace root (manifest discovery). */
    manifestCandidates(): string[];
    readText(relPath: string): Promise<string>;
    /** Write the manifest; returns how ('in-place' only under FSA). */
    saveManifest(relPath: string, text: string): Promise<'in-place' | 'download'>;
}

export function hasFsAccess(): boolean {
    return typeof window.showDirectoryPicker === 'function';
}

async function collectFsa(
    dir: FileSystemDirectoryHandle,
    prefix: string,
    out: { relPath: string; handle: FileSystemFileHandle }[],
): Promise<void> {
    for await (const [name, handle] of dir.entries()) {
        if (name.startsWith('.')) continue; // dot-entries are never cases
        const relPath = prefix ? `${prefix}/${name}` : name;
        if (handle.kind === 'file') {
            out.push({ relPath, handle });
        } else {
            await collectFsa(handle, relPath, out);
        }
    }
}

function fsaFile(relPath: string, handle: FileSystemFileHandle): WorkspaceFile {
    return {
        relPath,
        bytes: async () => (await handle.getFile()).arrayBuffer(),
        lastModified: async () => {
            const file = await handle.getFile();
            return Number.isFinite(file.lastModified) ? file.lastModified : null;
        },
    };
}

class FsaWorkspace implements Workspace {
    readonly kind = 'fsa' as const;
    readonly files: WorkspaceFile[];

    constructor(
        readonly root: FileSystemDirectoryHandle,
        entries: { relPath: string; handle: FileSystemFileHandle }[],
    ) {
        this.files = entries.map(({ relPath, handle }) => fsaFile(relPath, handle));
        this.byPath = new Map(entries.map((e) => [e.relPath, e.handle]));
    }

    private readonly byPath: Map<string, FileSystemFileHandle>;

    manifestCandidates(): string[] {
        return this.files
            .map((f) => f.relPath)
            .filter((p) => !p.includes('/') && p.endsWith('.json'));
    }

    async readText(relPath: string): Promise<string> {
        const handle = this.byPath.get(relPath);
        if (!handle) throw new Error(`meshio++: no workspace file '${relPath}'`);
        return (await handle.getFile()).text();
    }

    private async ensureWritable(): Promise<void> {
        const query = await this.root.queryPermission?.({ mode: 'readwrite' });
        if (query === 'granted') return;
        const granted = await this.root.requestPermission?.({ mode: 'readwrite' });
        if (granted !== 'granted') {
            throw new Error('meshio++: write permission to the directory was denied');
        }
    }

    async saveManifest(relPath: string, text: string): Promise<'in-place'> {
        await this.ensureWritable();
        // The manifest lives at the workspace root (it is the resolution
        // anchor for every relative source in it).
        const handle = await this.root.getFileHandle(relPath, { create: true });
        const writable = await handle.createWritable();
        await writable.write(text);
        await writable.close();
        if (!this.byPath.has(relPath)) {
            this.byPath.set(relPath, handle);
            this.files.push(fsaFile(relPath, handle));
        }
        return 'in-place';
    }
}

class FallbackWorkspace implements Workspace {
    readonly kind = 'fallback' as const;
    readonly files: WorkspaceFile[];

    constructor(fileList: File[], private readonly download: (text: string, name: string) => void) {
        // webkitRelativePath is `<picked-dir>/<...>`; strip the top segment so
        // paths are workspace-relative like the FSA side's.
        this.byPath = new Map(
            fileList
                .map((file) => {
                    const rel = file.webkitRelativePath || file.name;
                    const slash = rel.indexOf('/');
                    return [slash < 0 ? rel : rel.slice(slash + 1), file] as const;
                })
                .filter(([relPath]) => !relPath.split('/').some((s) => s.startsWith('.'))),
        );
        this.files = [...this.byPath.entries()].map(([relPath, file]) => ({
            relPath,
            bytes: () => file.arrayBuffer(),
            lastModified: async () =>
                Number.isFinite(file.lastModified) ? file.lastModified : null,
        }));
    }

    private readonly byPath: Map<string, File>;

    manifestCandidates(): string[] {
        return this.files
            .map((f) => f.relPath)
            .filter((p) => !p.includes('/') && p.endsWith('.json'));
    }

    async readText(relPath: string): Promise<string> {
        const file = this.byPath.get(relPath);
        if (!file) throw new Error(`meshio++: no workspace file '${relPath}'`);
        return file.text();
    }

    async saveManifest(relPath: string, text: string): Promise<'download'> {
        // The browser cannot write into the picked directory — a download is
        // the honest fallback, and the UI says so rather than pretending.
        this.download(text, relPath.split('/').pop() ?? 'dataset.json');
        return 'download';
    }
}

/**
 * Open a workspace from an existing directory handle — the picker's own
 * path and the persisted-handle "Reopen" restore both land here. A restored
 * handle arrives with permission `'prompt'`, so this MUST run under a user
 * gesture; denial is a named error the caller turns into clear-and-repick.
 */
export async function workspaceFromHandle(
    root: FileSystemDirectoryHandle,
): Promise<Workspace> {
    const query = (await root.queryPermission?.({ mode: 'readwrite' })) ?? 'granted';
    if (query !== 'granted') {
        const granted = await root.requestPermission?.({ mode: 'readwrite' });
        if (granted !== 'granted') {
            throw new Error(
                'meshio++: access to the saved directory was denied',
            );
        }
    }
    const entries: { relPath: string; handle: FileSystemFileHandle }[] = [];
    await collectFsa(root, '', entries);
    return new FsaWorkspace(root, entries);
}

/** Open a workspace with the File System Access API. */
export async function pickWorkspaceFsa(): Promise<Workspace> {
    if (!window.showDirectoryPicker) {
        throw new Error('meshio++: this browser has no File System Access API');
    }
    const root = await window.showDirectoryPicker({ mode: 'readwrite' });
    return workspaceFromHandle(root);
}

/** Open a workspace from a `webkitdirectory` input's FileList. */
export function workspaceFromFileList(
    fileList: FileList | File[],
    download: (text: string, name: string) => void,
): Workspace {
    return new FallbackWorkspace([...fileList], download);
}
