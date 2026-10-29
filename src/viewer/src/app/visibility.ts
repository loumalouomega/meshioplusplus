/**
 * Which surface cells are drawn.
 *
 * Several features hide cells independently (the selection tool's Hide and
 * Isolate, a hidden group). Each owns a named *layer* holding a mask of the
 * cells it leaves visible (1 = visible) or null when it hides nothing; a cell
 * is drawn only if every layer allows it. The composed mask is handed to the
 * renderer, and listeners (the selection overlay) are told when it changes.
 */
import type { Mask } from '../selection/mask';

export class Visibility {
    private readonly layers = new Map<string, Mask>();
    private readonly listeners = new Set<() => void>();
    private composed: Mask | null = null;
    private readonly apply: (mask: Mask | null) => void;

    constructor(apply: (mask: Mask | null) => void) {
        this.apply = apply;
    }

    /** Set a layer's mask, or pass null to clear it. */
    set(layer: string, mask: Mask | null): void {
        if (mask) this.layers.set(layer, mask);
        else this.layers.delete(layer);
        this.recompute();
    }

    get(layer: string): Mask | null {
        return this.layers.get(layer) ?? null;
    }

    /** Drop every layer (a new surface replaced the cells they described). */
    reset(): void {
        this.layers.clear();
        this.recompute();
    }

    /** Drop layers whose mask no longer matches the surface. */
    pruneTo(numCells: number): void {
        for (const [name, mask] of [...this.layers]) if (mask.length !== numCells) this.layers.delete(name);
        this.recompute();
    }

    /** The composed mask, or null when nothing is hidden. */
    combined(): Mask | null {
        return this.composed;
    }

    /** Number of cells currently hidden. */
    hiddenCount(): number {
        const m = this.composed;
        if (!m) return 0;
        let n = 0;
        for (let i = 0; i < m.length; i++) n += m[i] ^ 1;
        return n;
    }

    onChange(listener: () => void): () => void {
        this.listeners.add(listener);
        return () => this.listeners.delete(listener);
    }

    private recompute(): void {
        const masks = [...this.layers.values()];
        let next: Mask | null = null;
        if (masks.length === 1) next = masks[0];
        else if (masks.length > 1) {
            next = new Uint8Array(masks[0].length).fill(1);
            for (const m of masks) for (let i = 0; i < next.length; i++) next[i] &= m[i];
        }
        this.composed = next;
        this.apply(next);
        for (const listener of [...this.listeners]) listener();
    }
}
