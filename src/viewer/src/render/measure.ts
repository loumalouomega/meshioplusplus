/**
 * The measure tool's marks: the picked points and the segment between them.
 *
 * Plain vtk.js geometry in the scene, so it follows the camera and shows in
 * screenshots. The numbers are computed by {@link measureBetween}, which is
 * pure and unit-tested.
 */
import vtkActor from '@kitware/vtk.js/Rendering/Core/Actor';
import vtkCellArray from '@kitware/vtk.js/Common/Core/CellArray';
import vtkMapper from '@kitware/vtk.js/Rendering/Core/Mapper';
import vtkPolyData from '@kitware/vtk.js/Common/DataModel/PolyData';

import type { Vector3 } from '../types';
import type { Renderer } from './renderer';

export interface Measurement {
    distance: number;
    delta: Vector3;
}

/** Distance and per-axis difference from `a` to `b`. */
export function measureBetween(a: Vector3, b: Vector3): Measurement {
    const delta: Vector3 = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    return { distance: Math.hypot(...delta), delta };
}

const MARK_COLOR: Vector3 = [1.0, 0.35, 0.1];

export class MeasureOverlay {
    private readonly renderer: Renderer;
    private readonly mapper = vtkMapper.newInstance({ scalarVisibility: false });
    private readonly actor: vtkActor;

    constructor(renderer: Renderer) {
        this.renderer = renderer;
        this.actor = vtkActor.newInstance({ mapper: this.mapper });
        this.actor.setPickable(false);
        const prop = this.actor.getProperty();
        prop.setColor(...MARK_COLOR);
        prop.setLineWidth(3);
        prop.setPointSize(10);
        prop.setLighting(false);
        // Marks sit on the surface: pull them toward the viewer.
        this.mapper.setResolveCoincidentTopologyToPolygonOffset();
        this.mapper.setRelativeCoincidentTopologyLineOffsetParameters(-10, -10);
        this.mapper.setRelativeCoincidentTopologyPointOffsetParameters(-10, -10);
        this.actor.setVisibility(false);
        renderer.scene.addActor(this.actor);
    }

    /** Draw one or two points (and the segment when there are two). */
    set(points: readonly Vector3[]): void {
        if (points.length === 0) {
            this.clear();
            return;
        }
        const coords = new Float32Array(points.length * 3);
        points.forEach((p, i) => coords.set(p, i * 3));
        const poly = vtkPolyData.newInstance();
        poly.getPoints().setData(coords, 3);
        const verts = new Uint32Array(points.length * 2);
        points.forEach((_, i) => {
            verts[i * 2] = 1;
            verts[i * 2 + 1] = i;
        });
        poly.setVerts(vtkCellArray.newInstance({ values: verts }));
        if (points.length === 2) poly.setLines(vtkCellArray.newInstance({ values: Uint32Array.from([2, 0, 1]) }));
        this.mapper.setInputData(poly);
        this.actor.setVisibility(true);
        this.renderer.render();
    }

    clear(): void {
        this.actor.setVisibility(false);
        this.renderer.render();
    }
}
