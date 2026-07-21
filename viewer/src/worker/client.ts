/**
 * Main-thread half of the worker protocol.
 *
 * Only ever reached through a dynamic `import()` guarded by
 * `__VIEWER_EMBEDDED__`, so the wheel-bundled build drops it — and with it the
 * whole WASM binary, which that build has no use for.
 */
import type {
    ConvertResult,
    DistributiveOmit,
    InitResult,
    OpSpec,
    RenderResult,
    Request,
    Response,
} from './protocol';

type Pending = {
    resolve: (value: never) => void;
    reject: (reason: Error) => void;
    onProgress?: (stage: string) => void;
};

let worker: Worker | null = null;
let nextId = 1;
const pending = new Map<number, Pending>();

function ensureWorker(): Worker {
    if (worker) return worker;
    worker = new Worker(new URL('./worker.ts', import.meta.url), { type: 'module' });
    worker.onmessage = (event: MessageEvent<Response>) => {
        const message = event.data;
        const entry = pending.get(message.id);
        if (!entry) return;
        if (message.type === 'progress') {
            entry.onProgress?.(message.stage);
            return;
        }
        pending.delete(message.id);
        if (message.type === 'error') entry.reject(new Error(message.message));
        else entry.resolve(message as never);
    };
    worker.onerror = (e) => {
        // A worker-level failure (module load, WASM instantiation) never
        // resolves any single request, so fail them all rather than hang.
        const error = new Error(e.message || 'the meshio++ worker failed to start');
        for (const [, entry] of pending) entry.reject(error);
        pending.clear();
    };
    return worker;
}

function request<T>(
    message: DistributiveOmit<Request, 'id'>,
    options: { transfer?: Transferable[]; onProgress?: (stage: string) => void } = {}
): Promise<T> {
    const w = ensureWorker();
    const id = nextId++;
    return new Promise<T>((resolve, reject) => {
        pending.set(id, {
            resolve: resolve as (value: never) => void,
            reject,
            onProgress: options.onProgress,
        });
        w.postMessage({ id, ...message } as Request, options.transfer ?? []);
    });
}

/** Start the module and report which formats this build supports. */
export const init = (): Promise<InitResult> => request<InitResult>({ type: 'init' });

/** Stage `file` for the session and render it with an empty pipeline. */
export async function openFile(
    file: File,
    onProgress?: (stage: string) => void
): Promise<RenderResult> {
    const bytes = await file.arrayBuffer();
    return request<RenderResult>(
        { type: 'open', name: file.name, bytes },
        { transfer: [bytes], onProgress }
    );
}

/** Replay `ops` against the staged file. */
export const applyOps = (
    ops: OpSpec[],
    onProgress?: (stage: string) => void
): Promise<RenderResult> => request<RenderResult>({ type: 'apply', ops }, { onProgress });

/** Convert the staged file (with `ops` applied) to another format. */
export const convertTo = (outFormat: string, ops: OpSpec[]): Promise<ConvertResult> =>
    request<ConvertResult>({ type: 'convert', outFormat, ops });
