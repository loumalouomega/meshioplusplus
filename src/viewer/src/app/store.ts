/**
 * A tiny typed observable value.
 *
 * The viewer's features (selection, groups, theme, ...) each own one store and
 * publish changes through it, instead of main.ts reaching into every panel.
 * Dependency-free on purpose: it is bundled into unit tests unchanged.
 */
export type Listener<T> = (value: T, previous: T) => void;

export class Store<T> {
    private value: T;
    private readonly listeners = new Set<Listener<T>>();

    constructor(initial: T) {
        this.value = initial;
    }

    get(): T {
        return this.value;
    }

    /** Replace the value; listeners run only if it actually changed. */
    set(next: T): void {
        if (Object.is(next, this.value)) return;
        const previous = this.value;
        this.value = next;
        for (const listener of [...this.listeners]) listener(next, previous);
    }

    /** Update from the current value. */
    update(fn: (current: T) => T): void {
        this.set(fn(this.value));
    }

    /** Subscribe; returns the unsubscribe function. */
    subscribe(listener: Listener<T>): () => void {
        this.listeners.add(listener);
        return () => this.listeners.delete(listener);
    }
}
