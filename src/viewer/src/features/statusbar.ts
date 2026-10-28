/** The bottom status bar: what is open, how big it is, what is selected. */
import { maybe } from '../ui/dom';

export interface StatusInfo {
    file?: string;
    points?: number;
    cells?: number;
    selected?: number;
    hint?: string;
}

const fmt = (n: number) => n.toLocaleString('en-US');

function put(id: string, text: string): void {
    const el = maybe(id);
    if (el) el.textContent = text;
}

/** Update only the fields given; pass an empty string/0 to clear one. */
export function setStatusBar(info: StatusInfo): void {
    if (info.file !== undefined) put('sb-file', info.file);
    if (info.points !== undefined) put('sb-points', info.points ? `${fmt(info.points)} points` : '');
    if (info.cells !== undefined) put('sb-cells', info.cells ? `${fmt(info.cells)} surface cells` : '');
    if (info.selected !== undefined) put('sb-selection', info.selected ? `${fmt(info.selected)} selected` : '');
    if (info.hint !== undefined) put('sb-hint', info.hint);
}
