// meshio++ — MIT License (see LICENSE). Main authors: Vicente Mataix Ferrandiz
//
// vtk.js's own parse time (XML parse + base64 decode + inflate) of the VTP files
// `bench_vtp_inflate.py --keep DIR` writes (roadmap 3.4.5.1). It also records
// whether vtk.js reads each encoding at all: the uncompressed 4-byte-header file
// is expected to throw a RangeError for 8-byte types, which is why the viewer
// keeps the zlib default.
//
// Usage (vtk.js is a dependency of src/viewer; point --vtk at a directory whose
// node_modules holds @kitware/vtk.js@32.9.0, e.g. `npm install` there):
//
//   python benchmark/bench_vtp_inflate.py --keep /tmp/vtp
//   node benchmark/vtp_inflate.js --vtk src/viewer --dir /tmp/vtp [--repeats 5]

import { createRequire } from 'node:module';
import { readdirSync, readFileSync } from 'node:fs';
import { join, resolve } from 'node:path';
import { pathToFileURL } from 'node:url';

const args = process.argv.slice(2);
const opt = (name, fallback) => {
    const i = args.indexOf(`--${name}`);
    return i >= 0 ? args[i + 1] : fallback;
};
const vtkRoot = opt('vtk', 'src/viewer');
const dir = opt('dir');
const repeats = Number(opt('repeats', '5'));
if (!dir) {
    console.error('usage: node vtp_inflate.js --vtk <dir with node_modules> --dir <vtp dir>');
    process.exit(2);
}

const req = createRequire(join(resolve(vtkRoot), 'x.js'));
const readerUrl = pathToFileURL(
    req.resolve('@kitware/vtk.js/IO/XML/XMLPolyDataReader.js'),
).href;
const { default: XMLPolyDataReader } = await import(readerUrl);

const parse = (bytes) => {
    const reader = XMLPolyDataReader.newInstance();
    const buffer = bytes.buffer.slice(bytes.byteOffset, bytes.byteOffset + bytes.byteLength);
    reader.parseAsArrayBuffer(buffer);
    return reader.getOutputData(0);
};

console.log('file,bytes,status,min_ms,median_ms');
for (const name of readdirSync(dir).filter((f) => f.endsWith('.vtp')).sort()) {
    const bytes = readFileSync(join(dir, name));
    try {
        parse(bytes); // warm-up, and the read-at-all check
    } catch (e) {
        console.log(`${name},${bytes.length},${e.name}: ${e.message},,`);
        continue;
    }
    const times = [];
    for (let i = 0; i < repeats; i++) {
        const t0 = performance.now();
        parse(bytes);
        times.push(performance.now() - t0);
    }
    times.sort((a, b) => a - b);
    const median = times[Math.floor(times.length / 2)];
    console.log(`${name},${bytes.length},ok,${times[0].toFixed(2)},${median.toFixed(2)}`);
}
