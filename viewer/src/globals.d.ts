/// <reference types="vite/client" />
/** Ambient declarations: the test hook and Vite's compile-time constants. */
import type { ViewerState } from './types';

declare global {
    interface Window {
        /** See {@link ViewerState} — the documented Playwright hook. */
        __viewerState: ViewerState;
    }

    /**
     * True in the wheel-bundled build, where Python has already written the
     * mesh into the page and there is no WebAssembly. Vite `define`s it, so
     * Rollup eliminates the branches it guards.
     */
    const __VIEWER_EMBEDDED__: boolean;
}

export {};
