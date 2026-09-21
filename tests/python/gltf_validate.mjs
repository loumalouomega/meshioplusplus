// Runs the Khronos glTF-Validator (the `gltf-validator` npm package) over one
// file and prints its report as JSON: {"numErrors", "numWarnings", "numInfos",
// "messages": [...]}. Exit status 0 means no errors and no warnings.
//
//   node gltf_validate.mjs FILE.glb|FILE.gltf
//
// MESHIOPLUSPLUS_GLTF_VALIDATOR is a directory whose node_modules holds the
// package (e.g. the --prefix given to `npm install gltf-validator`); without it
// the package is resolved from the current directory's node_modules.
import { createRequire } from 'node:module';
import fs from 'node:fs';
import path from 'node:path';

const root = process.env.MESHIOPLUSPLUS_GLTF_VALIDATOR || process.cwd();
const require = createRequire(path.join(path.resolve(root), 'noop.js'));
const validator = require('gltf-validator');

const file = process.argv[2];
const dir = path.dirname(path.resolve(file));
const bytes = new Uint8Array(fs.readFileSync(file));

const report = await validator.validateBytes(bytes, {
    uri: path.basename(file),
    maxIssues: 0,
    externalResourceFunction: (uri) =>
        new Promise((resolve, reject) => {
            fs.readFile(path.join(dir, decodeURIComponent(uri)), (err, data) =>
                err ? reject(err) : resolve(new Uint8Array(data)));
        }),
});

const issues = report.issues;
console.log(JSON.stringify({
    numErrors: issues.numErrors,
    numWarnings: issues.numWarnings,
    numInfos: issues.numInfos,
    messages: issues.messages.map((m) => `${['error', 'warning', 'info', 'hint'][m.severity]} ${m.code} ${m.pointer}: ${m.message}`),
}));
process.exit(issues.numErrors === 0 && issues.numWarnings === 0 ? 0 : 1);
