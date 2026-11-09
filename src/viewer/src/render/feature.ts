/**
 * The measuring mode's highlight: the cells of the fitted features, drawn as
 * a translucent overlay in the scene (so it follows the camera and shows in
 * screenshots). The worker sends each feature's cells as plain triangles, so
 * no cell-id mapping between the fitted mesh and the drawn surface is needed.
 */
import vtkActor from '@kitware/vtk.js/Rendering/Core/Actor';
import vtkCellArray from '@kitware/vtk.js/Common/Core/CellArray';
import vtkMapper from '@kitware/vtk.js/Rendering/Core/Mapper';
import vtkPolyData from '@kitware/vtk.js/Common/DataModel/PolyData';

import type { Vector3 } from '../types';
import type { Renderer } from './renderer';

const FEATURE_COLOR: Vector3 = [0.1, 0.75, 0.45];
const SELECTED_COLOR: Vector3 = [1.0, 0.55, 0.1];

export class FeatureOverlay {
    private readonly renderer: Renderer;
    private readonly layers: { actor: vtkActor; mapper: vtkMapper }[] = [];

    constructor(renderer: Renderer) {
        this.renderer = renderer;
    }

    /** Draw each feature's triangles (nine floats each); `selected` ones in the accent colour. */
    set(features: readonly Float32Array[], selected: ReadonlySet<number> = new Set()): void {
        while (this.layers.length < features.length) {
            const mapper = vtkMapper.newInstance({ scalarVisibility: false });
            mapper.setResolveCoincidentTopologyToPolygonOffset();
            mapper.setRelativeCoincidentTopologyPolygonOffsetParameters(-2, -2);
            const actor = vtkActor.newInstance({ mapper });
            actor.setPickable(false);
            actor.getProperty().setOpacity(0.55);
            this.renderer.scene.addActor(actor);
            this.layers.push({ actor, mapper });
        }
        this.layers.forEach((layer, i) => {
            const tris = features[i];
            if (!tris || tris.length < 9) {
                layer.actor.setVisibility(false);
                return;
            }
            const n = tris.length / 9;
            const poly = vtkPolyData.newInstance();
            poly.getPoints().setData(tris, 3);
            const cells = new Uint32Array(n * 4);
            for (let t = 0; t < n; ++t) cells.set([3, 3 * t, 3 * t + 1, 3 * t + 2], 4 * t);
            poly.setPolys(vtkCellArray.newInstance({ values: cells }));
            layer.mapper.setInputData(poly);
            layer.actor.getProperty().setColor(...(selected.has(i) ? SELECTED_COLOR : FEATURE_COLOR));
            layer.actor.setVisibility(true);
        });
        this.renderer.render();
    }

    clear(): void {
        this.set([]);
    }
}
