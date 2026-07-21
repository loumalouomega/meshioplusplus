/**
 * Mapping data arrays onto colours.
 *
 * Owns the colour-by vocabulary: which arrays are colourable, the stable key
 * that identifies one choice, and the mapper/LUT configuration for it.
 *
 * The point-vs-cell branch in {@link configureScalarMapper} is the pair of
 * settings most viewers get wrong: point data interpolates before mapping
 * (smooth gradients across a face), cell data must not (flat per-cell colour).
 */
import type vtkDataArray from '@kitware/vtk.js/Common/Core/DataArray';
import type vtkMapper from '@kitware/vtk.js/Rendering/Core/Mapper';
import type vtkPolyData from '@kitware/vtk.js/Common/DataModel/PolyData';
import type vtkColorTransferFunction from '@kitware/vtk.js/Rendering/Core/ColorTransferFunction';
import { ColorMode, ScalarMode } from '@kitware/vtk.js/Rendering/Core/Mapper/Constants';
import { VectorMode } from '@kitware/vtk.js/Common/Core/ScalarsToColors/Constants';

import type { ArrayEntry, ScalarRange } from '../types';

/** The "no data, just show the shape" entry of the colour-by menu. */
export const SOLID_COLOR = '';

/**
 * Arrays whose names start with this are plumbing, not something to colour by.
 * `surface:parent_cell` is kept in the polydata for the picker's provenance
 * row but must never reach the menu.
 */
const INTERNAL_PREFIX = 'surface:';

/** Stable identifier for one colour-by choice. */
export function colorKey(entry: ArrayEntry): string {
    return `${entry.location}:${entry.name}:${entry.component}`;
}

/**
 * Every colourable array, one entry per renderable component.
 *
 * A multi-component array yields a magnitude entry plus one per component:
 * exactly what the lossless `convertSurface` path preserves, and what a
 * `readMesh`/`writeMesh` round-trip would have flattened away.
 */
export function listArrays(polydata: vtkPolyData | null): ArrayEntry[] {
    if (!polydata) return [];
    const out: ArrayEntry[] = [];
    const sources = [
        ['point', polydata.getPointData()],
        ['cell', polydata.getCellData()],
    ] as const;

    for (const [location, attrs] of sources) {
        for (let i = 0; i < attrs.getNumberOfArrays(); i++) {
            const array = attrs.getArrayByIndex(i);
            if (!array) continue;
            const name = array.getName();
            if (name.startsWith(INTERNAL_PREFIX)) continue;
            const n = array.getNumberOfComponents();
            if (n === 1) {
                out.push({ location, name, component: -1, label: name });
                continue;
            }
            out.push({ location, name, component: -1, label: `${name} (magnitude)` });
            for (let c = 0; c < n; c++) {
                out.push({ location, name, component: c, label: `${name}[${c}]` });
            }
        }
    }
    return out;
}

/**
 * The value range of one component, and how many entries were non-finite.
 *
 * vtk.js's own `getRange` already skips NaN, but it seeds `min`/`max` from
 * `±Number.MAX_VALUE` and returns those untouched when *every* value is
 * non-finite — an inverted range that maps the whole mesh to one end of the
 * colormap. Compute it here so that case is detectable and reportable.
 */
export function computeRange(array: vtkDataArray, component: number): ScalarRange {
    const values = array.getData();
    const stride = array.getNumberOfComponents();
    const magnitude = component < 0 && stride > 1;
    const offset = component < 0 ? 0 : component;

    let min = Number.POSITIVE_INFINITY;
    let max = Number.NEGATIVE_INFINITY;
    let nanCount = 0;
    const tuples = values.length / stride;

    for (let t = 0; t < tuples; t++) {
        let v: number;
        if (magnitude) {
            let sum = 0;
            for (let c = 0; c < stride; c++) {
                const x = values[t * stride + c];
                sum += x * x;
            }
            v = Math.sqrt(sum);
        } else {
            v = values[t * stride + offset];
        }
        if (!Number.isFinite(v)) {
            nanCount++;
            continue;
        }
        if (v < min) min = v;
        if (v > max) max = v;
    }

    if (min > max) return { min: 0, max: 0, nanCount };
    return { min, max, nanCount };
}

/**
 * Point a mapper at one array.
 *
 * ## Why `interpolateScalarsBeforeMapping` is off
 *
 * Textbook advice is to turn it on for point data — map each vertex to a
 * colour, then interpolate colours — because interpolating the *scalar* first
 * can walk through colours the data never had. In vtk.js 32.9.0 it is
 * measurably wrong for a range spanning zero: on a field over [-0.149, 0.149]
 * with a diverging map it renders **entirely in the warm half** (0 blue pixels
 * out of 275k coloured ones), and reordering the range/mapper setup does not
 * help. Presumably its scalar-to-texture-coordinate path mishandles a negative
 * range.
 *
 * So it stays off for both associations. The visible cost is small — colours
 * are still interpolated across a face, just after mapping instead of before —
 * and the alternative is silently showing half the colormap.
 *
 * For cell data it must be off regardless: interpolating a per-cell value
 * across its own cell would smear a categorical tag into a meaningless
 * gradient.
 */
export function configureScalarMapper(
    mapper: vtkMapper,
    lut: vtkColorTransferFunction,
    entry: ArrayEntry
): void {
    mapper.setScalarVisibility(true);
    mapper.setColorMode(ColorMode.MAP_SCALARS);
    mapper.setColorByArrayName(entry.name);
    if (entry.location === 'point') {
        mapper.setScalarMode(ScalarMode.USE_POINT_FIELD_DATA);
        mapper.setInterpolateScalarsBeforeMapping(false);
    } else {
        mapper.setScalarMode(ScalarMode.USE_CELL_FIELD_DATA);
        mapper.setInterpolateScalarsBeforeMapping(false);
    }
    lut.setVectorComponent(entry.component);
    lut.setVectorMode(
        entry.component < 0 ? VectorMode.MAGNITUDE : VectorMode.COMPONENT
    );
}

/** The array behind an entry, or null if it is gone (e.g. after an operation). */
export function arrayOf(
    polydata: vtkPolyData | null,
    entry: ArrayEntry
): vtkDataArray | null {
    if (!polydata) return null;
    const attrs =
        entry.location === 'point' ? polydata.getPointData() : polydata.getCellData();
    return attrs.getArrayByName(entry.name);
}
