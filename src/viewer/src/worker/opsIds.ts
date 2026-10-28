/**
 * Which pipeline operations keep the surface's cells (same cells, same order).
 *
 * Selections and groups are masks over surface cells, so they stay valid only
 * through operations that neither add, remove nor reorder cells. This is a
 * small pure module shared by the client and, later, the worker, so both agree.
 *
 * Measured against the bundled block mesh (`quality`, `smooth`, `partition`
 * leave `surface:parent_cell` identical). `clean` is treated as destructive
 * even though it can be a no-op, because it may drop degenerate or duplicate
 * cells; `refine`, `section`, `isosurface`, `voxelize` and `computeSdf` build
 * new cells.
 */
import type { OpSpec } from './protocol';

const PRESERVING: ReadonlySet<string> = new Set(['quality', 'smooth', 'partition', 'gradient']);

export function opPreservesCellIds(op: OpSpec): boolean {
    return PRESERVING.has(op.op);
}

export function opsPreserveCellIds(ops: readonly OpSpec[]): boolean {
    return ops.every(opPreservesCellIds);
}
