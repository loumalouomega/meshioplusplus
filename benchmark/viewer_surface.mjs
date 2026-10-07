// meshio++ — MIT License (see LICENSE). Main authors: Vicente Mataix Ferrandiz
//
// The browser viewer's render path, old against new (roadmap 3.4.5), on the
// same WASM build:
//
//   vtp     convertSurfaceOps -> /out.vtp (zlib) -> FS.readFile, then vtk.js's
//           XMLPolyDataReader on the main thread (XML, base64, inflate);
//   buffers surfaceBuffersOps -> typed arrays, then a vtkPolyData built from
//           them as Renderer.loadSurface builds it.
//
// Each is split into its worker half and its main-thread half, and `read` is a
// plain readMesh of the same file, an upper bound on the read every apply
// repeats (it also copies the mesh into JS). Before timing, the two polydata
// are compared point for point, cell for cell and array for array, so a row
// is only printed for a file both paths render identically.
//
// Usage (vtk.js is a dependency of src/viewer; `npm install` there first):
//
//   python benchmark/bench_vtp_inflate.py --keep /tmp/vs   # surfaces, any mesh files do
//   node benchmark/viewer_surface.mjs --vtk src/viewer --dir /tmp/vs [--repeats 5] [--variant seq|mt]

import { createRequire } from 'node:module';
import { readdirSync, readFileSync } from 'node:fs';
import { join, resolve } from 'node:path';
import { pathToFileURL } from 'node:url';

import { loadMeshioPlusPlus } from '../src/wasm/src/index.mjs';

const args = process.argv.slice(2);
const opt = (name, fallback) => {
    const i = args.indexOf(`--${name}`);
    return i >= 0 ? args[i + 1] : fallback;
};
const vtkRoot = opt('vtk', 'src/viewer');
const dir = opt('dir');
const repeats = Number(opt('repeats', '5'));
const variant = opt('variant', 'seq');
const only = opt('only', '');
if (!dir) {
    console.error('usage: node viewer_surface.mjs --vtk <dir with node_modules> --dir <mesh dir>');
    process.exit(2);
}

const req = createRequire(join(resolve(vtkRoot), 'x.js'));
const vtk = async (path) => (await import(pathToFileURL(req.resolve(`@kitware/vtk.js/${path}`)).href)).default;
const XMLPolyDataReader = await vtk('IO/XML/XMLPolyDataReader.js');
const vtkPolyData = await vtk('Common/DataModel/PolyData.js');
const vtkCellArray = await vtk('Common/Core/CellArray.js');
const vtkDataArray = await vtk('Common/Core/DataArray.js');

const m = await loadMeshioPlusPlus({}, { variant });

const median = (xs) => [...xs].sort((a, b) => a - b)[Math.floor(xs.length / 2)];
function time(fn) {
    fn(); // warm-up
    const ts = [];
    for (let i = 0; i < repeats; i++) {
        const t0 = performance.now();
        fn();
        ts.push(performance.now() - t0);
    }
    return median(ts);
}

// --- the two paths, each in its worker and main-thread halves ------------- //

const vtpWorker = (path) => {
    m.convertSurfaceOps(path, '/out.vtp', [], { keepProvenance: true });
    const bytes = m.FS.readFile('/out.vtp', { encoding: 'binary' });
    m.FS.unlink('/out.vtp');
    return bytes.buffer;
};
const vtpMain = (buffer) => {
    const reader = XMLPolyDataReader.newInstance();
    reader.parseAsArrayBuffer(buffer);
    return reader.getOutputData(0);
};
const buffersWorker = (path) => m.surfaceBuffersOps(path, [], { keepProvenance: true }).surface;
const buffersMain = (s) => {
    // Renderer.loadSurface, verbatim but for the class wrapper.
    const polydata = vtkPolyData.newInstance();
    polydata.getPoints().setData(s.points, 3);
    if (s.verts.length) polydata.setVerts(vtkCellArray.newInstance({ values: s.verts }));
    if (s.lines.length) polydata.setLines(vtkCellArray.newInstance({ values: s.lines }));
    if (s.polys.length) polydata.setPolys(vtkCellArray.newInstance({ values: s.polys }));
    for (const a of s.pointData)
        polydata.getPointData().addArray(
            vtkDataArray.newInstance({ name: a.name, values: a.values, numberOfComponents: a.components }),
        );
    for (const a of s.cellData)
        polydata.getCellData().addArray(
            vtkDataArray.newInstance({ name: a.name, values: a.values, numberOfComponents: a.components }),
        );
    return polydata;
};

// --- equivalence ---------------------------------------------------------- //

const same = (a, b, what) => {
    if (a.length !== b.length) throw new Error(`${what}: length ${a.length} != ${b.length}`);
    for (let i = 0; i < a.length; i++)
        if (Number(a[i]) !== Number(b[i])) throw new Error(`${what}[${i}]: ${a[i]} != ${b[i]}`);
};
function assertSamePolyData(p, q) {
    same(p.getPoints().getData(), q.getPoints().getData(), 'points');
    for (const sec of ['Verts', 'Lines', 'Polys'])
        same(p[`get${sec}`]().getData(), q[`get${sec}`]().getData(), sec);
    for (const loc of ['getPointData', 'getCellData']) {
        const a = p[loc](),
            b = q[loc]();
        if (a.getNumberOfArrays() !== b.getNumberOfArrays()) throw new Error(`${loc}: array count`);
        for (let i = 0; i < a.getNumberOfArrays(); i++) {
            const x = a.getArrayByIndex(i),
                y = b.getArrayByIndex(i);
            if (x.getName() !== y.getName()) throw new Error(`${loc}: ${x.getName()} != ${y.getName()}`);
            if (x.getNumberOfComponents() !== y.getNumberOfComponents())
                throw new Error(`${x.getName()}: components`);
            same(x.getData(), y.getData(), x.getName());
        }
    }
}

// --- run ------------------------------------------------------------------ //

try {
    m.FS.mkdir('/in');
} catch {
    // already there
}
console.log(
    'file,points,cells,vtp_bytes,read_ms,vtp_worker_ms,vtp_main_ms,vtp_total_ms,buffers_worker_ms,buffers_main_ms,buffers_total_ms',
);
for (const name of readdirSync(dir).sort()) {
    if (only && !name.includes(only)) continue;
    const path = `/in/${name}`;
    m.FS.writeFile(path, readFileSync(join(dir, name)));
    let reference;
    try {
        reference = vtpMain(vtpWorker(path));
    } catch (e) {
        console.error(`skip ${name}: ${e.message}`);
        m.FS.unlink(path);
        continue;
    }
    assertSamePolyData(reference, buffersMain(buffersWorker(path)));

    const vtpBytes = vtpWorker(path).byteLength;
    const read = time(() => m.readMesh(path));
    const vw = time(() => vtpWorker(path));
    const vtpBuffer = vtpWorker(path);
    const vm = time(() => vtpMain(vtpBuffer));
    const bw = time(() => buffersWorker(path));
    const surface = buffersWorker(path);
    const bm = time(() => buffersMain(surface));
    const f = (x) => x.toFixed(1);
    console.log(
        [
            name,
            reference.getNumberOfPoints(),
            reference.getNumberOfCells(),
            vtpBytes,
            f(read),
            f(vw),
            f(vm),
            f(vw + vm),
            f(bw),
            f(bm),
            f(bw + bm),
        ].join(','),
    );
    m.FS.unlink(path);
}
