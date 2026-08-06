/**
 * Ambient declarations for the dataset page: the `__datasetState` test hook
 * and the File System Access API surface TypeScript's DOM lib does not yet
 * ship (`showDirectoryPicker` and the async directory iterator live behind
 * `DOM.AsyncIterable`, which this tsconfig's `lib` does not include).
 * Declared here — type-checked declaration merging, never `as any`.
 */

import type { DatasetState } from './types';

declare global {
    interface Window {
        __datasetState: DatasetState;
        showDirectoryPicker?: (options?: {
            mode?: 'read' | 'readwrite';
        }) => Promise<FileSystemDirectoryHandle>;
    }

    interface FileSystemDirectoryHandle {
        entries(): AsyncIterableIterator<
            [string, FileSystemDirectoryHandle | FileSystemFileHandle]
        >;
        queryPermission?(descriptor: {
            mode: 'read' | 'readwrite';
        }): Promise<PermissionState>;
        requestPermission?(descriptor: {
            mode: 'read' | 'readwrite';
        }): Promise<PermissionState>;
    }
}

export {};
