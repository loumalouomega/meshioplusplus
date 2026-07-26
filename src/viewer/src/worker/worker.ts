/**
 * The meshio++ WASM worker.
 *
 * Every read, operation and conversion happens here so the main thread stays
 * responsive on a large file. Payloads move as transferable `ArrayBuffer`s, so
 * nothing is copied across the boundary.
 *
 * The session model is the important part. The **original file bytes stay
 * staged in MEMFS** for as long as the file is open, and every pipeline change
 * replays the whole operation list against them through one `convertSurfaceOps`
 * call. That is what makes undo exact — pop an operation, replay, done — with
 * no inverse operations and no snapshots. It also means `apply` must *not*
 * unlink the staged file, which is the one place this deviates from the
 * "always clean up after yourself" discipline everywhere else in here.
 */
import { loadMeshioPlusPlus } from '@meshioplusplus/wasm';
// Vite emits the binary as a hashed asset and hands back a URL already
// prefixed with `base`. Doing this explicitly, rather than letting the
// Emscripten glue resolve its own .wasm, is what makes the module load
// correctly from inside a module worker *and* through the copied
// `file:../wasm` dependency.
//
// Both variants are imported: the loader auto-selects the threaded build
// (meshioplusplus_wasm_mt) when this page is cross-origin isolated, else the
// sequential one, and each `.mjs` requests the matching `.wasm` by name --
// which locateFile below resolves. Vite tree-shakes neither, so both hashed
// assets ship; the browser only fetches the one that is actually loaded.
import wasmUrl from '@meshioplusplus/wasm/dist/meshioplusplus_wasm.wasm?url';
import wasmMtUrl from '@meshioplusplus/wasm/dist/meshioplusplus_wasm_mt.wasm?url';

import type { Vector3 } from '../types';
import type {
    MeshMeta,
    OpReport,
    OpSpec,
    Request,
    Response,
} from './protocol';

type Module = Awaited<ReturnType<typeof loadMeshioPlusPlus>>;

let mio: Module | null = null;

const IN_DIR = '/in';
const SURFACE_PATH = '/out.vtp';

interface Session {
    path: string;
    name: string;
    format: string;
    meta: MeshMeta;
    bounds: { min: Vector3; max: Vector3 } | null;
}

let session: Session | null = null;

async function ready(): Promise<Module> {
    if (!mio) {
        // locateFile receives the requested filename; return the URL matching
        // whichever variant the loader picked (the _mt glue asks for the _mt
        // binary). Passing a fixed URL would hand the threaded glue the
        // sequential binary and fail to instantiate.
        mio = await loadMeshioPlusPlus({
            locateFile: (path: string) => (path.includes('_mt') ? wasmMtUrl : wasmUrl),
        });
        try {
            mio.FS.mkdir(IN_DIR);
        } catch {
            // already there
        }
    }
    return mio;
}

function unlink(m: Module, path: string): void {
    try {
        m.FS.unlink(path);
    } catch {
        // MEMFS never reclaims on its own, so this matters after several
        // files -- but a missing file is not worth failing a request over.
    }
}

/** Read a file out of the virtual FS as a detached, transferable buffer. */
function take(m: Module, path: string): ArrayBuffer {
    // Emscripten's readFile decodes to a string unless told otherwise, which
    // would corrupt any binary format on the way out.
    const data = m.FS.readFile(path, { encoding: 'binary' });
    if (typeof data === 'string') {
        throw new Error(`meshio++: ${path} was read as text, not bytes`);
    }
    // .slice() so we never transfer a view onto HEAPU8: that would detach the
    // WASM heap and kill the module for every later request.
    return data.slice().buffer as ArrayBuffer;
}

function boundsOf(meta: MeshMeta): { min: Vector3; max: Vector3 } | null {
    const { bboxMin, bboxMax } = meta;
    if (!bboxMin || !bboxMax || bboxMin.length < 3 || bboxMax.length < 3) return null;
    return {
        min: [bboxMin[0], bboxMin[1], bboxMin[2]],
        max: [bboxMax[0], bboxMax[1], bboxMax[2]],
    };
}

/** Replay `ops` against the staged file and return the renderable surface. */
function renderPipeline(m: Module, ops: OpSpec[]): { vtp: ArrayBuffer; report: OpReport } {
    if (!session) throw new Error('no mesh is open');
    const report = m.convertSurfaceOps(session.path, SURFACE_PATH, ops, {
        inFormat: session.format,
        // The picker needs the provenance array; the colour-by menu filters it
        // out by name, so keeping it costs nothing visible.
        keepProvenance: true,
    }) as OpReport;
    const vtp = take(m, SURFACE_PATH);
    unlink(m, SURFACE_PATH);
    return { vtp, report };
}

function post(message: Response, transfer: Transferable[] = []): void {
    (self as unknown as Worker).postMessage(message, transfer);
}

async function handle(request: Request): Promise<void> {
    const m = await ready();
    const { id } = request;

    switch (request.type) {
        case 'init': {
            const formats =
                typeof m.availableFormats === 'function'
                    ? m.availableFormats()
                    : { readers: [], writers: [] };
            post({ id, type: 'result', kind: 'init', formats, backend: m.meshBackend() });
            return;
        }

        case 'open': {
            // One staged file at a time; the previous session's bytes go now.
            if (session) unlink(m, session.path);
            session = null;

            // The extension is how the format is inferred, so the original
            // name has to survive -- but a path separator in it would escape
            // the staging directory.
            const safe = request.name.replace(/[/\\]/g, '_') || 'mesh';
            const path = `${IN_DIR}/${safe}`;
            m.FS.writeFile(path, new Uint8Array(request.bytes));

            let format = '';
            try {
                format = m.sniffFormat(path) || '';
            } catch {
                format = '';
            }

            post({ id, type: 'progress', stage: 'read' });
            const meta = m.readMetadata(path, format) as unknown as MeshMeta;
            session = {
                path,
                name: request.name,
                format: meta.format || format,
                meta,
                bounds: boundsOf(meta),
            };

            post({ id, type: 'progress', stage: 'surface' });
            const { vtp, report } = renderPipeline(m, []);
            post(
                {
                    id,
                    type: 'result',
                    kind: 'render',
                    vtp,
                    meta,
                    format: session.format,
                    report,
                    bounds: session.bounds,
                },
                [vtp]
            );
            return;
        }

        case 'apply': {
            if (!session) throw new Error('no mesh is open');
            post({ id, type: 'progress', stage: 'operations' });
            // Deliberately no unlink: the staged bytes are the pipeline's
            // input every time, which is what makes undo exact.
            const { vtp, report } = renderPipeline(m, request.ops);
            post(
                {
                    id,
                    type: 'result',
                    kind: 'render',
                    vtp,
                    meta: session.meta,
                    format: session.format,
                    report,
                    bounds: session.bounds,
                },
                [vtp]
            );
            return;
        }

        case 'convert': {
            if (!session) throw new Error('no mesh is open');
            const outPath = `/converted.${request.outFormat}`;
            try {
                if (request.ops.length === 0) {
                    // No operations: convert the file as it is, which keeps
                    // volume cells rather than exporting only the surface.
                    m.convert(session.path, outPath, {
                        inFormat: session.format,
                        outFormat: request.outFormat,
                    });
                } else {
                    m.convertSurfaceOps(session.path, outPath, request.ops, {
                        inFormat: session.format,
                        outFormat: request.outFormat,
                    });
                }
                const bytes = take(m, outPath);
                post(
                    {
                        id,
                        type: 'result',
                        kind: 'convert',
                        bytes,
                        outFormat: request.outFormat,
                    },
                    [bytes]
                );
            } finally {
                unlink(m, outPath);
            }
            return;
        }

        case 'close': {
            if (session) unlink(m, session.path);
            session = null;
            post({ id, type: 'result', kind: 'closed' });
            return;
        }
    }
}

self.onmessage = async (event: MessageEvent<Request>) => {
    try {
        await handle(event.data);
    } catch (e) {
        // -fwasm-exceptions plus the bindings' EM_ASM rethrow means meshio++
        // read/write failures arrive here as real Errors with real messages.
        const message = e instanceof Error ? e.message : String(e);
        post({ id: event.data.id, type: 'error', message });
    }
};
