/**
 * A small SVG line chart for training curves — hand-rolled rather than a
 * charting dependency: the viewer has exactly two runtime dependencies
 * (vtk.js and the meshio++ WASM package) and a loss curve is two polylines
 * and two axes.
 *
 * Conventions follow this repo's data-visualization rules: ONE y-axis
 * (never a second scale), series colours assigned in fixed order from the
 * categorical slots (validated against this page's own dark surface —
 * blue/orange, worst-pair CVD ΔE 26.8), a legend AND direct end-labels so
 * identity is never colour alone, recessive hairline grid, 2px marks, and
 * a hover crosshair with a tooltip (a chart in a browser is interactive).
 * Axis text wears the muted ink token, never a series colour.
 *
 * The tick helpers are pure and unit-tested (`tests/viewer/unit/chart.test.mjs`).
 */

/** Categorical slots 1-2 for the dark surface — train, then valid. */
export const SERIES_COLORS = ['#3987e5', '#d95926'];

export interface ChartSeries {
    name: string;
    color: string;
    /** `[x, y]` pairs, ascending in x; non-finite y values are skipped. */
    points: [number, number][];
}

export interface ChartOptions {
    /** Log-scale the y axis (losses span orders of magnitude). */
    logY?: boolean;
    xLabel?: string;
    yLabel?: string;
    /** Element the hover tooltip is written into (positioned by the caller's
     * CSS; hidden when the pointer leaves the plot). */
    tooltip?: HTMLElement | null;
    /** Formats a hovered x value in the tooltip (default: the raw number). */
    formatX?: (x: number) => string;
}

const SVG_NS = 'http://www.w3.org/2000/svg';
const WIDTH = 640;
const HEIGHT = 300;
const MARGIN = { left: 62, right: 74, top: 12, bottom: 34 };

/** Round numbers spanning [min, max]: 1/2/5 x 10^k, at most `count` of them. */
export function niceTicks(min: number, max: number, count = 5): number[] {
    if (!Number.isFinite(min) || !Number.isFinite(max)) return [];
    if (min === max) return [min];
    const span = max - min;
    const rough = span / Math.max(count, 1);
    const magnitude = 10 ** Math.floor(Math.log10(rough));
    const step = [1, 2, 5, 10].map((m) => m * magnitude).find((s) => s >= rough) ?? magnitude * 10;
    const ticks: number[] = [];
    for (let v = Math.ceil(min / step) * step; v <= max + step / 1e6; v += step) {
        ticks.push(Number(v.toPrecision(12)));
    }
    return ticks;
}

/** Decade ticks over [min, max] (positive), thinned to at most `count`. */
export function logTicks(min: number, max: number, count = 6): number[] {
    if (!(min > 0) || !(max > 0) || min > max) return [];
    const lo = Math.floor(Math.log10(min));
    const hi = Math.ceil(Math.log10(max));
    const decades: number[] = [];
    // `10 ** -4` is 0.00009999999999999999 in IEEE double; round the decade
    // back to its clean value so a tick reads as the number it names.
    for (let e = lo; e <= hi; e += 1) decades.push(Number((10 ** e).toPrecision(15)));
    const inside = decades.filter((v) => v >= min / 10 && v <= max * 10);
    const stride = Math.max(1, Math.ceil(inside.length / Math.max(count, 1)));
    return inside.filter((_, i) => i % stride === 0);
}

/** Compact tick text: 1.2e-3, 0.25, 12, 3.4k. */
export function formatTick(value: number): string {
    if (!Number.isFinite(value)) return '';
    const magnitude = Math.abs(value);
    if (magnitude === 0) return '0';
    if (magnitude < 1e-3 || magnitude >= 1e5) return value.toExponential(1).replace('e+', 'e');
    return String(Number(value.toPrecision(3)));
}

function el<K extends keyof SVGElementTagNameMap>(
    tag: K,
    attrs: Record<string, string | number> = {},
): SVGElementTagNameMap[K] {
    const node = document.createElementNS(SVG_NS, tag);
    for (const [key, value] of Object.entries(attrs)) node.setAttribute(key, String(value));
    return node;
}

interface Scale {
    (value: number): number;
    invert(pixel: number): number;
}

function linearScale(d0: number, d1: number, r0: number, r1: number): Scale {
    const span = d1 - d0 || 1;
    const scale = ((value: number) => r0 + ((value - d0) / span) * (r1 - r0)) as Scale;
    scale.invert = (pixel: number) => d0 + ((pixel - r0) / (r1 - r0 || 1)) * span;
    return scale;
}

function logScale(d0: number, d1: number, r0: number, r1: number): Scale {
    const l0 = Math.log10(d0);
    const l1 = Math.log10(d1);
    const inner = linearScale(l0, l1 || l0 + 1, r0, r1);
    const scale = ((value: number) => inner(Math.log10(Math.max(value, Number.MIN_VALUE)))) as Scale;
    scale.invert = (pixel: number) => 10 ** inner.invert(pixel);
    return scale;
}

/** Draw `series` into `svg`, replacing whatever was there. */
export function renderLineChart(
    svg: SVGSVGElement,
    series: ChartSeries[],
    options: ChartOptions = {},
): void {
    svg.replaceChildren();
    svg.setAttribute('viewBox', `0 0 ${WIDTH} ${HEIGHT}`);
    svg.setAttribute('preserveAspectRatio', 'none');
    const drawable = series
        .map((s) => ({
            ...s,
            points: s.points.filter(([x, y]) => Number.isFinite(x) && Number.isFinite(y)),
        }))
        .filter((s) => s.points.length > 0);
    if (!drawable.length) return;

    const xs = drawable.flatMap((s) => s.points.map(([x]) => x));
    const ys = drawable.flatMap((s) => s.points.map(([, y]) => y));
    const positive = ys.filter((y) => y > 0);
    const logY = !!options.logY && positive.length === ys.length && positive.length > 0;
    const xMin = Math.min(...xs);
    const xMax = Math.max(...xs);
    let yMin = Math.min(...ys);
    let yMax = Math.max(...ys);
    if (yMin === yMax) {
        // A flat series still needs a band to draw in.
        yMin = logY ? yMin / 2 : yMin - 1;
        yMax = logY ? yMax * 2 : yMax + 1;
    }
    const plotLeft = MARGIN.left;
    const plotRight = WIDTH - MARGIN.right;
    const plotTop = MARGIN.top;
    const plotBottom = HEIGHT - MARGIN.bottom;
    const x = linearScale(xMin, xMax === xMin ? xMin + 1 : xMax, plotLeft, plotRight);
    const y = logY
        ? logScale(yMin, yMax, plotBottom, plotTop)
        : linearScale(yMin, yMax, plotBottom, plotTop);

    // Grid + axes: hairlines, recessive, behind the marks.
    const yTicks = logY ? logTicks(yMin, yMax) : niceTicks(yMin, yMax);
    for (const tick of yTicks) {
        const py = y(tick);
        if (py < plotTop - 0.5 || py > plotBottom + 0.5) continue;
        svg.append(
            el('line', {
                x1: plotLeft,
                x2: plotRight,
                y1: py,
                y2: py,
                stroke: 'var(--line)',
                'stroke-width': 1,
            }),
        );
        const label = el('text', {
            x: plotLeft - 8,
            y: py + 4,
            'text-anchor': 'end',
            fill: 'var(--muted)',
            'font-size': 11,
        });
        label.textContent = formatTick(tick);
        svg.append(label);
    }
    const xTicks = niceTicks(xMin, xMax, 6).filter((t) => Number.isInteger(t));
    for (const tick of xTicks.length ? xTicks : [xMin, xMax]) {
        const px = x(tick);
        const label = el('text', {
            x: px,
            y: plotBottom + 18,
            'text-anchor': 'middle',
            fill: 'var(--muted)',
            'font-size': 11,
        });
        label.textContent = formatTick(tick);
        svg.append(label);
    }
    svg.append(
        el('line', {
            x1: plotLeft,
            x2: plotRight,
            y1: plotBottom,
            y2: plotBottom,
            stroke: 'var(--muted)',
            'stroke-width': 1,
        }),
    );
    if (options.xLabel) {
        const label = el('text', {
            x: (plotLeft + plotRight) / 2,
            y: HEIGHT - 4,
            'text-anchor': 'middle',
            fill: 'var(--muted)',
            'font-size': 11,
        });
        label.textContent = options.xLabel;
        svg.append(label);
    }
    if (options.yLabel) {
        const label = el('text', {
            x: 12,
            y: (plotTop + plotBottom) / 2,
            'text-anchor': 'middle',
            fill: 'var(--muted)',
            'font-size': 11,
            transform: `rotate(-90 12 ${(plotTop + plotBottom) / 2})`,
        });
        label.textContent = options.yLabel;
        svg.append(label);
    }

    // Marks, then a direct end-label per series (identity never colour alone).
    for (const s of drawable) {
        const path = s.points.map(([px, py]) => `${x(px)},${y(py)}`).join(' ');
        svg.append(
            el('polyline', {
                points: path,
                fill: 'none',
                stroke: s.color,
                'stroke-width': 2,
                'stroke-linejoin': 'round',
                'stroke-linecap': 'round',
            }),
        );
        const last = s.points[s.points.length - 1]!;
        const label = el('text', {
            x: Math.min(x(last[0]) + 8, WIDTH - 4),
            y: y(last[1]) + 4,
            fill: s.color,
            'font-size': 11,
        });
        label.textContent = s.name;
        svg.append(label);
    }

    // Hover: a crosshair rule plus a per-series dot, and the tooltip text.
    const crosshair = el('line', {
        y1: plotTop,
        y2: plotBottom,
        stroke: 'var(--muted)',
        'stroke-width': 1,
        'stroke-dasharray': '3 3',
        visibility: 'hidden',
    });
    svg.append(crosshair);
    const dots = drawable.map((s) => {
        const dot = el('circle', { r: 4, fill: s.color, visibility: 'hidden' });
        svg.append(dot);
        return dot;
    });
    const hit = el('rect', {
        x: plotLeft,
        y: plotTop,
        width: Math.max(plotRight - plotLeft, 1),
        height: Math.max(plotBottom - plotTop, 1),
        fill: 'transparent',
    });
    svg.append(hit);

    const tooltip = options.tooltip ?? null;
    const formatX = options.formatX ?? ((value: number) => formatTick(value));
    const hide = () => {
        crosshair.setAttribute('visibility', 'hidden');
        for (const dot of dots) dot.setAttribute('visibility', 'hidden');
        if (tooltip) tooltip.hidden = true;
    };
    hit.addEventListener('pointerleave', hide);
    hit.addEventListener('pointermove', (event) => {
        const box = svg.getBoundingClientRect();
        if (!box.width || !box.height) return;
        // The viewBox scales; map client pixels back into user units.
        const userX = ((event.clientX - box.left) / box.width) * WIDTH;
        const target = x.invert(Math.min(Math.max(userX, plotLeft), plotRight));
        crosshair.setAttribute('x1', String(x(target)));
        crosshair.setAttribute('x2', String(x(target)));
        crosshair.setAttribute('visibility', 'visible');
        const rows: string[] = [];
        drawable.forEach((s, i) => {
            let nearest = s.points[0]!;
            for (const point of s.points) {
                if (Math.abs(point[0] - target) < Math.abs(nearest[0] - target)) nearest = point;
            }
            const dot = dots[i]!;
            dot.setAttribute('cx', String(x(nearest[0])));
            dot.setAttribute('cy', String(y(nearest[1])));
            dot.setAttribute('visibility', 'visible');
            rows.push(`${s.name} ${formatTick(nearest[1])}`);
        });
        if (tooltip) {
            const nearestX = drawable[0]!.points.reduce((best, point) =>
                Math.abs(point[0] - target) < Math.abs(best[0] - target) ? point : best,
            );
            tooltip.textContent = `${formatX(nearestX[0])} · ${rows.join(' · ')}`;
            tooltip.hidden = false;
        }
    });
}
