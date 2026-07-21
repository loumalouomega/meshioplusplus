/**
 * Click-to-inspect: what is this cell, and what are its values?
 *
 * `vtkCellPicker` rather than `vtkHardwareSelector`: the latter has no `.d.ts`,
 * needs an extra render pass plus a `readPixels` GPU stall, and is fragile on
 * the software renderer CI runs on. `CellPicker` is a CPU ray/cell
 * intersection over the polydata — a few milliseconds for a boundary surface
 * of a few hundred thousand triangles, which is fine for something that fires
 * once per click.
 *
 * Two guards instead of an acceleration structure, which would be premature
 * for a once-per-click interaction: pick **only on click, never on move**, and
 * refuse above {@link MAX_PICKABLE_CELLS}.
 */
import vtkActor from '@kitware/vtk.js/Rendering/Core/Actor';
import vtkCellArray from '@kitware/vtk.js/Common/Core/CellArray';
import vtkCellPicker from '@kitware/vtk.js/Rendering/Core/CellPicker';
import vtkDataArray from '@kitware/vtk.js/Common/Core/DataArray';
import vtkMapper from '@kitware/vtk.js/Rendering/Core/Mapper';
import vtkPolyData from '@kitware/vtk.js/Common/DataModel/PolyData';

import type { PickInfo, PickValue, Vector3 } from '../types';
import { Representation, type Renderer } from './renderer';

/** Above this, inspecting is disabled rather than made slow. */
export const MAX_PICKABLE_CELLS = 2_000_000;

/** The provenance array `convertSurfaceOps(keepProvenance: true)` leaves behind. */
const PARENT_ARRAY = 'surface:parent_cell';

/** High-contrast against every shipped colormap. */
const HIGHLIGHT_COLOR: Vector3 = [1.0, 0.85, 0.1];

export class Picker {
    private readonly renderer: Renderer;
    private readonly picker = vtkCellPicker.newInstance();
    private readonly highlightMapper = vtkMapper.newInstance({ scalarVisibility: false });
    private readonly highlightActor: vtkActor;

    constructor(renderer: Renderer) {
        this.renderer = renderer;
        this.highlightActor = vtkActor.newInstance({ mapper: this.highlightMapper });
        const prop = this.highlightActor.getProperty();
        prop.setRepresentation(Representation.WIREFRAME);
        prop.setColor(...HIGHLIGHT_COLOR);
        prop.setLineWidth(3);
        // The highlight lies exactly on the surface it outlines; without an
        // offset toward the viewer it z-fights into a dashed mess.
        this.highlightMapper.setResolveCoincidentTopologyToPolygonOffset();
        this.highlightMapper.setRelativeCoincidentTopologyLineOffsetParameters(-8, -8);
        this.highlightActor.setVisibility(false);
        renderer.scene.addActor(this.highlightActor);
    }

    /** Whether this mesh is small enough to inspect. */
    canPick(): boolean {
        const data = this.renderer.data;
        return !!data && data.getNumberOfCells() <= MAX_PICKABLE_CELLS;
    }

    clear(): void {
        this.highlightActor.setVisibility(false);
        this.renderer.render();
    }

    /**
     * Pick at a canvas position.
     *
     * @param x canvas-relative x, in CSS pixels
     * @param y canvas-relative y, in CSS pixels (browser convention, y-down)
     * @param height canvas height in CSS pixels, for the y flip
     */
    pickAt(x: number, y: number, height: number): PickInfo | null {
        const polydata = this.renderer.data;
        if (!polydata || !this.canPick()) return null;

        // Display space is y-up; the browser gives y-down.
        this.picker.pick([x, height - y, 0], this.renderer.scene);
        if (this.picker.getActors().length === 0) {
            this.clear();
            return null;
        }

        const cellId = this.picker.getCellId();
        if (cellId < 0 || cellId >= polydata.getNumberOfCells()) {
            this.clear();
            return null;
        }

        const position = this.picker.getMapperPosition();
        const pointIds = cellPointIds(polydata, cellId);
        const pointId = nearestPoint(polydata, pointIds, position);

        this.highlight(polydata, pointIds);

        return {
            cellId,
            cellType: cellTypeName(polydata, cellId),
            pointId,
            position: pointCoordinates(polydata, pointId),
            parentCell: parentCellOf(polydata, cellId),
            cellValues: valuesAt(polydata.getCellData(), cellId),
            pointValues: valuesAt(polydata.getPointData(), pointId),
        };
    }

    /** Outline just the picked cell, as its own tiny polydata. */
    private highlight(polydata: vtkPolyData, pointIds: number[]): void {
        const source = polydata.getPoints().getData();
        const points = new Float64Array(pointIds.length * 3);
        pointIds.forEach((id, i) => {
            points[i * 3] = source[id * 3];
            points[i * 3 + 1] = source[id * 3 + 1];
            points[i * 3 + 2] = source[id * 3 + 2];
        });

        const cell = new Uint32Array(pointIds.length + 1);
        cell[0] = pointIds.length;
        for (let i = 0; i < pointIds.length; i++) cell[i + 1] = i;

        const outline = vtkPolyData.newInstance();
        outline.getPoints().setData(points, 3);
        outline.setPolys(vtkCellArray.newInstance({ values: cell }));

        this.highlightMapper.setInputData(outline);
        this.highlightActor.setVisibility(true);
        this.renderer.render();
    }
}

// --- polydata helpers ----------------------------------------------------- //

function cellPointIds(polydata: vtkPolyData, cellId: number): number[] {
    const info = polydata.getCellPoints(cellId) as {
        cellPointIds?: ArrayLike<number>;
    } | null;
    return info?.cellPointIds ? Array.from(info.cellPointIds) : [];
}

/**
 * The cell's vertex closest to where the ray actually hit.
 *
 * `CellPicker` declares no point getter, so rather than reach for an
 * undeclared one this computes it from the connectivity — deterministic, and
 * it cannot break when vtk.js changes what it does or does not expose.
 */
function nearestPoint(
    polydata: vtkPolyData,
    pointIds: number[],
    position: Vector3
): number {
    const coords = polydata.getPoints().getData();
    let best = pointIds[0] ?? -1;
    let bestDistance = Infinity;
    for (const id of pointIds) {
        const dx = coords[id * 3] - position[0];
        const dy = coords[id * 3 + 1] - position[1];
        const dz = coords[id * 3 + 2] - position[2];
        const d = dx * dx + dy * dy + dz * dz;
        if (d < bestDistance) {
            bestDistance = d;
            best = id;
        }
    }
    return best;
}

function pointCoordinates(polydata: vtkPolyData, pointId: number): Vector3 {
    if (pointId < 0) return [0, 0, 0];
    const c = polydata.getPoints().getData();
    return [c[pointId * 3], c[pointId * 3 + 1], c[pointId * 3 + 2]];
}

/** meshio++ cell-type name for a polydata cell, by its vertex count. */
function cellTypeName(polydata: vtkPolyData, cellId: number): string {
    const n = cellPointIds(polydata, cellId).length;
    if (n === 3) return 'triangle';
    if (n === 4) return 'quad';
    if (n === 2) return 'line';
    if (n === 1) return 'vertex';
    return `polygon${n}`;
}

/**
 * The original cell this facet came from, if that is still meaningful.
 *
 * Only when the operation pipeline is empty: `surface:parent_cell` names a cell
 * of whatever mesh the boundary was extracted from, so after a refine or a
 * section it refers to the *post-operation* mesh rather than the file. Chaining
 * the ids through every operation is possible but genuinely hard for marginal
 * value, so the UI says which it means instead of pretending.
 */
function parentCellOf(polydata: vtkPolyData, cellId: number): number | null {
    const array = polydata.getCellData().getArrayByName(PARENT_ARRAY);
    if (!array) return null;
    const value = array.getData()[cellId];
    return Number.isFinite(value) ? Number(value) : null;
}

/** Every array's value at one element, all components. */
function valuesAt(
    attributes: {
        getNumberOfArrays(): number;
        // vtk.js types this Nullable, and it really can be null for a sparse
        // attribute slot -- hence the guard in the loop.
        getArrayByIndex(i: number): vtkDataArray | null;
    },
    index: number
): PickValue[] {
    if (index < 0) return [];
    const out: PickValue[] = [];
    for (let i = 0; i < attributes.getNumberOfArrays(); i++) {
        const array = attributes.getArrayByIndex(i);
        if (!array) continue;
        const name = array.getName();
        // Provenance is reported on its own row, not as a data array.
        if (name === PARENT_ARRAY) continue;
        const n = array.getNumberOfComponents();
        const data = array.getData();
        const components: number[] = [];
        for (let c = 0; c < n; c++) components.push(Number(data[index * n + c]));
        out.push({ name, components });
    }
    return out;
}
