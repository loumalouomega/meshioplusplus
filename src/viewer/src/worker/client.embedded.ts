/**
 * Build-time stub replacing `client.ts` in the embedded (wheel) build.
 *
 * That build renders a VTP that Python already produced, so it needs no WASM
 * at all. Aliasing the module out is not merely a size optimization: Vite's
 * worker plugin emits the `new Worker(new URL('./worker.ts', ...))` chunk when
 * it *transforms* the real client, whether or not the code is reachable, so
 * dead-code elimination alone still shipped a multi-megabyte worker beside a
 * page that could never call it.
 *
 * Nothing here should ever run — `main.ts` returns before reaching it under
 * `__VIEWER_EMBEDDED__` — so these throw rather than fail quietly.
 */
import type { ConvertResult, InitResult, RenderResult } from './protocol';

function unavailable(): never {
    throw new Error(
        'this build of the meshio++ viewer renders an embedded mesh only; ' +
            'use the online viewer to open other files'
    );
}

export const init = (): Promise<InitResult> => unavailable();
export const openFile = (): Promise<RenderResult> => unavailable();
export const applyOps = (): Promise<RenderResult> => unavailable();
export const convertTo = (): Promise<ConvertResult> => unavailable();
