import { createRequire } from 'node:module';
import { fileURLToPath } from 'node:url';

import { defineConfig } from 'vite';
import { viteSingleFile } from 'vite-plugin-singlefile';

/** Wrapper functions the viewer calls; see src/wasm/src/index.mjs. */
const REQUIRED_WASM_API = [
    'readMesh',
    'readMetadata',
    'convert',
    'convertSurface',
    'convertSurfaceOps',
    'sniffFormat',
    'availableFormats',
    'meshBackend',
];

const REFRESH_INSTRUCTIONS =
    'Rebuild the package and refresh the copy:\n' +
    '  source /path/to/emsdk/emsdk_env.sh\n' +
    '  ./build/configure-wasm.sh --build        # from the repository root\n' +
    '  cd viewer && rm -rf node_modules/@meshioplusplus && npm install\n' +
    "(plain `npm install` is not enough -- npm sees no dependency change and skips the copy)";

/**
 * Fail the build if the installed meshio++ package is not the one in src/wasm/.
 *
 * npm *copies* a `file:` dependency rather than symlinking it, and worse, a
 * later `npm install` will not refresh that copy: nothing in package.json
 * changed, so npm skips it. A viewer built after a `src/wasm/` rebuild therefore
 * bundles the *previous* module and misbehaves in ways that look like viewer
 * bugs -- a missing binding surfaces as "m.convertSurface is not a function",
 * and a changed binding surfaces as nothing at all.
 *
 * Comparing the wrapper's exports catches the first; comparing the binary
 * catches the second, which is the one that cost real debugging time.
 */
function checkWasmPackage() {
    return {
        name: 'meshioplusplus:check-wasm-package',
        async buildStart() {
            const fs = await import('node:fs/promises');
            const require = createRequire(import.meta.url);

            let wrapper;
            try {
                wrapper = require.resolve('@meshioplusplus/wasm/src/index.mjs');
            } catch {
                this.error(
                    '@meshioplusplus/wasm is not installed. Run `npm install` in viewer/.'
                );
            }

            const source = await fs.readFile(wrapper, 'utf8');
            const missing = REQUIRED_WASM_API.filter((name) => !source.includes(`${name}:`));
            if (missing.length) {
                this.error(
                    `the installed @meshioplusplus/wasm does not forward: ` +
                        `${missing.join(', ')}.\n${REFRESH_INSTRUCTIONS}`
                );
            }

            // Both artifacts must be present and current: the sequential
            // meshioplusplus_wasm and the threaded meshioplusplus_wasm_mt (the
            // loader auto-selects between them, and the worker imports both).
            for (const binary of [
                'dist/meshioplusplus_wasm.wasm',
                'dist/meshioplusplus_wasm_mt.wasm',
            ]) {
                const installed = fileURLToPath(
                    new URL(`./node_modules/@meshioplusplus/wasm/${binary}`, import.meta.url)
                );
                const built = fileURLToPath(new URL(`../wasm/${binary}`, import.meta.url));
                const [a, b] = await Promise.all([
                    fs.readFile(installed).catch(() => null),
                    fs.readFile(built).catch(() => null),
                ]);
                if (!b) {
                    this.error(
                        `${built} does not exist -- src/wasm/dist/ is gitignored, so the package ` +
                            `must be built before the viewer (build both variants; do not pass ` +
                            `--seq-only).\n${REFRESH_INSTRUCTIONS}`
                    );
                }
                if (!a || !a.equals(b)) {
                    this.error(
                        `the installed @meshioplusplus/wasm ${binary} differs from the one in ` +
                            `src/wasm/dist/.\n${REFRESH_INSTRUCTIONS}`
                    );
                }
            }
        },
    };
}

// Two builds from one source tree.
//
//   build:web    the public demo, deployed under /meshioplusplus/viewer/. Reads
//                any format via the meshio++ WASM package in a worker.
//   build:embed  a single self-contained HTML file shipped inside the Python
//                wheel and used by `meshioplusplus.view(backend="browser")`.
//                There, Python has the real C++ core and has already written
//                the VTP, so the browser only has to render: the WASM worker is
//                behind a `__VIEWER_EMBEDDED__` guard and Rollup drops it,
//                along with the ~2 MB .wasm binary. Being a single file is also
//                what makes it work over file:// with no local HTTP server.
const embedded = !!process.env.VITE_VIEWER_EMBEDDED;

/**
 * Inject the COOP/COEP service worker (public/coi-serviceworker.js) into the
 * web build's index.html only. It makes GitHub Pages cross-origin isolated so
 * the worker can load the THREADED wasm artifact; without it the loader
 * auto-falls-back to the sequential one, so this is a speed enhancement, never
 * a correctness requirement. It must NOT reach the embed build: that build is a
 * single self-contained file with no wasm and no `public/`, and an external
 * <script src> would break the one-file guarantee. The tag is injected here
 * rather than written into index.html precisely so the embed build never sees
 * it. The relative src resolves against the document URL (i.e. under `base`),
 * so the service worker's scope covers the whole viewer on Pages.
 */
function injectCoiServiceWorker() {
    return {
        name: 'meshioplusplus:coi-serviceworker',
        transformIndexHtml() {
            if (embedded) return;
            return [
                {
                    tag: 'script',
                    attrs: { src: 'coi-serviceworker.js' },
                    injectTo: 'head-prepend',
                },
            ];
        },
    };
}

export default defineConfig({
    base: process.env.VITE_BASE ?? '/meshioplusplus/viewer/',
    plugins: embedded
        ? [viteSingleFile({ removeViteModuleLoader: true })]
        : [checkWasmPackage(), injectCoiServiceWorker()],
    define: { __VIEWER_EMBEDDED__: JSON.stringify(embedded) },
    resolve: {
        alias: embedded
            ? {
                  // Replace the WASM client with a stub so the worker chunk is
                  // never emitted. Dead-code elimination is not enough here:
                  // Vite's worker plugin emits the `new Worker(new URL(...))`
                  // chunk when it *transforms* the module, reachable or not,
                  // which shipped a 5.6 MB worker beside a page that could
                  // never call it.
                  // The key must match the import specifier in main.ts
                  // *exactly* -- Vite aliases are literal string matches, and
                  // a trailing `.ts` here silently matches nothing, which puts
                  // the whole 5.7 MB worker back into the single-file build.
                  './worker/client': fileURLToPath(
                      new URL('./src/worker/client.embedded.ts', import.meta.url)
                  ),
              }
            : {},
    },
    // The embedded build must be exactly one file; the samples belong to the
    // public demo only.
    publicDir: embedded ? false : 'public',
    // The meshio++ package is a `file:../wasm` dependency, so npm symlinks it
    // and the dev server has to be allowed to read outside the project root.
    server: { fs: { allow: ['..'] } },
    worker: { format: 'es' },
    build: {
        target: 'es2022',
        // The embed build must inline every asset; nothing may be fetched.
        assetsInlineLimit: embedded ? 100_000_000 : 4096,
        chunkSizeWarningLimit: 4096,
    },
});
