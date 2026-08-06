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
    ArraySummary,
    MeshMeta,
    OpReport,
    OpSpec,
    PlanEntry,
    Request,
    Response,
    StagedSource,
    StageResult,
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

/**
 * The dataset page's staged entry — a second, independent slot beside the
 * viewer's `session`; the two never interact. Exactly ONE entry is resident
 * at a time (MEMFS never reclaims, so staging a manifest's every case would
 * exhaust the 32-bit heap): `stageEntry` evicts the previous one first.
 */
interface DatasetStage {
    dir: string;
    plan: PlanEntry[];
    /** Per-plan-step inferred format ('' = let the reader infer). */
    formats: string[];
}

let datasetStage: DatasetStage | null = null;
let stageCounter = 0;

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

// --- dataset staging (dataset.html) --------------------------------------- //

/** Create every directory of `path` (a `mkdir -p` over the typed FS API). */
function mkdirs(m: Module, path: string): void {
    const parts = path.split('/').filter(Boolean);
    let current = '';
    for (const part of parts) {
        current += `/${part}`;
        try {
            m.FS.mkdir(current);
        } catch {
            // already there
        }
    }
}

// The wasm package's FS type enumerates only the four calls the viewer uses;
// everything else sits behind its `[key: string]: unknown` index signature.
// Narrow, named casts (the vtk-shims convention), never `any`.
function fsReaddir(m: Module, dir: string): string[] {
    return (m.FS.readdir as (path: string) => string[])(dir);
}

function fsRmdir(m: Module, dir: string): void {
    (m.FS.rmdir as (path: string) => void)(dir);
}

/** Remove a staging directory recursively (MEMFS never reclaims on its own). */
function removeTree(m: Module, dir: string): void {
    let names: string[];
    try {
        names = fsReaddir(m, dir).filter((n) => n !== '.' && n !== '..');
    } catch {
        return; // not there
    }
    for (const name of names) {
        const path = `${dir}/${name}`;
        try {
            m.FS.unlink(path); // a file
        } catch {
            removeTree(m, path); // a directory
        }
    }
    try {
        fsRmdir(m, dir);
    } catch {
        // best effort
    }
}

function sniffOf(m: Module, path: string): string {
    try {
        return m.sniffFormat(path) || '';
    } catch {
        return '';
    }
}

function stageEntry(
    m: Module,
    files: { relPath: string; bytes: ArrayBuffer }[],
    source: StagedSource,
): StageResult {
    if (datasetStage) removeTree(m, datasetStage.dir);
    datasetStage = null;

    stageCounter += 1;
    const dir = `/data/e${stageCounter}`;
    for (const file of files) {
        const segments = file.relPath.split('/');
        if (
            file.relPath.startsWith('/') ||
            segments.some((s) => s === '..' || s === '')
        ) {
            throw new Error(`meshio++: unsafe staged path '${file.relPath}'`);
        }
        // Directory structure is PRESERVED here (unlike `open`'s `/`→`_`
        // flattening): a manifest Pattern's directory part is literal, so the
        // staged tree must mirror the workspace's.
        mkdirs(m, `${dir}/${segments.slice(0, -1).join('/')}`);
        m.FS.writeFile(`${dir}/${file.relPath}`, new Uint8Array(file.bytes));
    }

    // The wasm side is the single authority on plan ordering and times —
    // identical to Python's `DatasetEntry.entries()`.
    const prefix = (p: string) => `${dir}/${p}`;
    const memfsSource: string | string[] = source.Paths
        ? source.Paths.map(prefix)
        : prefix(source.Pattern ?? source.Path ?? '');
    let plan = m.sequenceEntries(memfsSource, {
        format: source.Format,
        times: source.Times,
        timeFrom: source.TimeFrom,
        sort: source.Sort,
    }) as PlanEntry[];

    // A multi-step file's steps cannot be previewed directly — the render
    // path (`convertSurfaceOps`) has no step selector — so fan it out once
    // into per-step files; every later scrub tick is then the same cheap
    // single-file render every other Source kind uses.
    if (plan.some((entry) => entry.step !== 0)) {
        mkdirs(m, `${dir}/steps`);
        const inPath = plan[0].path;
        const fanned = m.timeseriesToSequence(
            inPath,
            `${dir}/steps/s_{step}.vtu`,
            source.Format,
        );
        plan = plan.map((entry, i) => ({ ...entry, path: fanned[i] ?? entry.path, step: 0 }));
    }

    const formats = plan.map((entry) => source.Format ?? sniffOf(m, entry.path));
    const meta = m.readMetadata(plan[0].path, formats[0]) as unknown as MeshMeta;
    datasetStage = { dir, plan, formats };
    return { plan, meta };
}

function stagedStep(step: number): { path: string; format: string } {
    if (!datasetStage) throw new Error('no dataset entry is staged');
    const entry = datasetStage.plan[step];
    if (!entry) throw new Error(`step ${step} is out of range`);
    return { path: entry.path, format: datasetStage.formats[step] ?? '' };
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

        case 'stageEntry': {
            post({ id, type: 'progress', stage: 'stage' });
            const result = stageEntry(m, request.files, request.source);
            post({ id, type: 'result', kind: 'staged', ...result });
            return;
        }

        case 'previewStep': {
            const { path, format } = stagedStep(request.step);
            post({ id, type: 'progress', stage: 'surface' });
            const meta = m.readMetadata(path, format) as unknown as MeshMeta;
            const report = m.convertSurfaceOps(path, SURFACE_PATH, [], {
                inFormat: format,
                keepProvenance: true,
            }) as OpReport;
            const vtp = take(m, SURFACE_PATH);
            unlink(m, SURFACE_PATH);
            post(
                {
                    id,
                    type: 'result',
                    kind: 'render',
                    vtp,
                    meta,
                    format: meta.format || format,
                    report,
                    bounds: boundsOf(meta),
                },
                [vtp],
            );
            return;
        }

        case 'summaryStep': {
            const { path, format } = stagedStep(request.step);
            post({ id, type: 'progress', stage: 'summary' });
            // Mesh-taking is fine for scalar aggregates (never a rendering
            // path — the flat JS Mesh cannot carry multi-component arrays,
            // but their component stats survive via the sidecar convention).
            let mesh = m.readMeshSelective(path, { format });
            try {
                // quality:* rows ride the same ArraySummary shape. Their NaN
                // means "metric N/A for this cell type" BY DESIGN — the UI's
                // bad-case rules must exclude quality-prefixed names.
                mesh = m.attachQuality(mesh);
            } catch {
                // quality is best-effort; the plain summary still answers
            }
            const arrays = m.dataInfo(mesh).map(
                (info) =>
                    ({
                        location: info.location,
                        name: info.name,
                        dtype: info.dtype,
                        numValues: info.numValues,
                        numComponents: info.numComponents,
                        min: info.min,
                        max: info.max,
                        mean: info.mean,
                        numNan: info.numNan,
                        numInf: info.numInf,
                    }) satisfies ArraySummary,
            );
            post({ id, type: 'result', kind: 'summary', arrays });
            return;
        }

        case 'evictEntry': {
            if (datasetStage) removeTree(m, datasetStage.dir);
            datasetStage = null;
            post({ id, type: 'result', kind: 'evicted' });
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
