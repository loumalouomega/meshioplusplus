/**
 * Selection highlight: one flat-coloured actor per selection set, built from
 * the cells in the set's mask (points shared with the surface).
 *
 * A polygon offset pulls the overlay toward the viewer so it does not z-fight
 * with the surface it sits on.
 */
import vtkActor from '@kitware/vtk.js/Rendering/Core/Actor';
import vtkCellArray from '@kitware/vtk.js/Common/Core/CellArray';
import vtkMapper from '@kitware/vtk.js/Rendering/Core/Mapper';
import vtkPolyData from '@kitware/vtk.js/Common/DataModel/PolyData';

import type { Renderer } from '../render/renderer';
import { subsetCells } from '../render/subset';
import type { Vector3 } from '../types';
import type { Mask } from './mask';

interface Entry {
    actor: vtkActor;
    mapper: vtkMapper;
}

export class SelectionOverlay {
    private readonly renderer: Renderer;
    private readonly entries = new Map<number, Entry>();

    constructor(renderer: Renderer) {
        this.renderer = renderer;
    }

    /** Draw `mask` for set `id` in `color`; an empty or null mask removes it. */
    set(id: number, color: Vector3, mask: Mask | null, visible: Mask | null): void {
        const data = this.renderer.data;
        const geometry = this.renderer.surfaceGeometry();
        if (!data || !geometry || !mask) {
            this.remove(id);
            return;
        }

        // Hidden cells are not drawn, so neither is their highlight.
        let effective = mask;
        if (visible) {
            effective = new Uint8Array(mask.length);
            for (let i = 0; i < mask.length; i++) effective[i] = mask[i] & visible[i];
        }

        let any = false;
        for (let i = 0; i < effective.length && !any; i++) any = effective[i] === 1;
        if (!any) {
            this.remove(id);
            return;
        }

        const get = (c: { getData(): ArrayLike<number> } | null | undefined) =>
            Uint32Array.from(c ? c.getData() : []);
        const subset = subsetCells(
            get(data.getVerts()),
            get(data.getLines()),
            get(data.getPolys()),
            effective
        );
        const poly = vtkPolyData.newInstance();
        poly.setPoints(data.getPoints());
        if (subset.verts.length) poly.setVerts(vtkCellArray.newInstance({ values: subset.verts }));
        if (subset.lines.length) poly.setLines(vtkCellArray.newInstance({ values: subset.lines }));
        if (subset.polys.length) poly.setPolys(vtkCellArray.newInstance({ values: subset.polys }));

        let entry = this.entries.get(id);
        if (!entry) {
            const mapper = vtkMapper.newInstance({ scalarVisibility: false });
            mapper.setResolveCoincidentTopologyToPolygonOffset();
            mapper.setRelativeCoincidentTopologyPolygonOffsetParameters(-3, -3);
            mapper.setRelativeCoincidentTopologyLineOffsetParameters(-6, -6);
            this.renderer.makeClippable(mapper);
            const actor = vtkActor.newInstance({ mapper });
            actor.setPickable(false);
            const prop = actor.getProperty();
            // Lit, so the set keeps the surface's shading; a strong ambient term
            // keeps its own colour readable at any angle.
            prop.setAmbient(0.55);
            prop.setDiffuse(0.55);
            prop.setLineWidth(3);
            prop.setPointSize(6);
            this.renderer.scene.addActor(actor);
            entry = { actor, mapper };
            this.entries.set(id, entry);
        }
        entry.actor.getProperty().setColor(...color);
        entry.mapper.setInputData(poly);
        entry.actor.setVisibility(true);
        this.renderer.render();
    }

    remove(id: number): void {
        const entry = this.entries.get(id);
        if (!entry) return;
        this.renderer.scene.removeActor(entry.actor);
        this.renderer.forgetClippable(entry.mapper);
        entry.actor.delete();
        this.entries.delete(id);
        this.renderer.render();
    }

    clear(): void {
        for (const id of [...this.entries.keys()]) this.remove(id);
    }
}
