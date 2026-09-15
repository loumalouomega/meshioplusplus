/**
 * Best-effort persistence of the picked directory handle across reloads.
 *
 * A real `FileSystemDirectoryHandle` is structured-cloneable, so IndexedDB
 * can store it; on the next visit the page offers a "Reopen" button whose
 * click (a user gesture — `requestPermission` demands one) re-grants access
 * without re-picking.
 *
 * **Every function here is best-effort internally.** A store that cannot
 * work — no `indexedDB`, private mode, quota, or a `DataCloneError` because
 * the "handle" is a plain mock object (the e2e FSA mock's handles carry
 * function properties and cannot clone) — degrades to exactly today's
 * behaviour: `saveHandle` resolves anyway, `loadHandle` resolves `null`,
 * `clearHandle` resolves. Callers never need a try/catch.
 *
 * `window.indexedDB` is read at call time, which is also the test seam: the
 * Playwright spec swaps in a Map-backed fake via `addInitScript`.
 */

const DB_NAME = 'meshioplusplus-dataset';
/** One object store for everything (the directory handle, the scan cache,
 * the server connection): a second store would need a DB version bump and a
 * migration, for no benefit — keys are namespaced instead. */
const STORE = 'handles';
const KEY = 'workspace-root';

function withStore<T>(
    mode: IDBTransactionMode,
    fn: (store: IDBObjectStore) => IDBRequest<T>,
): Promise<T | null> {
    return new Promise((resolve) => {
        let request: IDBOpenDBRequest;
        try {
            request = window.indexedDB.open(DB_NAME, 1);
        } catch {
            resolve(null);
            return;
        }
        request.onupgradeneeded = () => {
            request.result.createObjectStore(STORE);
        };
        request.onerror = () => resolve(null);
        request.onsuccess = () => {
            const db = request.result;
            try {
                const tx = db.transaction(STORE, mode);
                const op = fn(tx.objectStore(STORE));
                op.onsuccess = () => {
                    resolve((op.result as T) ?? null);
                    db.close();
                };
                op.onerror = () => {
                    resolve(null);
                    db.close();
                };
            } catch {
                // e.g. DataCloneError from a non-cloneable mock handle
                resolve(null);
                db.close();
            }
        };
    });
}

/** Remember the workspace root; failures are silently absorbed. */
export async function saveHandle(handle: FileSystemDirectoryHandle): Promise<void> {
    await withStore('readwrite', (store) => store.put(handle, KEY));
}

/** The previously saved root, or null (absent, unsupported, or failed). */
export async function loadHandle(): Promise<FileSystemDirectoryHandle | null> {
    const value = await withStore<unknown>('readonly', (store) => store.get(KEY));
    if (
        value &&
        typeof value === 'object' &&
        (value as FileSystemDirectoryHandle).kind === 'directory'
    ) {
        return value as FileSystemDirectoryHandle;
    }
    return null;
}

/** Forget the saved root (a denied/stale handle must not keep reappearing). */
export async function clearHandle(): Promise<void> {
    await withStore('readwrite', (store) => store.delete(KEY));
}

/** Store any structured-cloneable value under a namespaced key
 * (`scan:...`, `server`); failures are silently absorbed. */
export async function putValue(key: string, value: unknown): Promise<void> {
    await withStore('readwrite', (store) => store.put(value, key));
}

/** Read a value stored with {@link putValue}, or null. */
export async function getValue<T>(key: string): Promise<T | null> {
    const value = await withStore<unknown>('readonly', (store) => store.get(key));
    return value === undefined ? null : (value as T | null);
}
