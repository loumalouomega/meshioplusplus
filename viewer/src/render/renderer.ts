/**
 * The vtk.js scene.
 *
 * vtk.js's only mesh data model is `vtkPolyData` — it ships no unstructured
 * grid and no XML unstructured-grid reader — which is why everything reaching
 * this module is VTP, and why a volume mesh is shown by its boundary surface.
 * The meshio++ WASM side does that extraction, so nothing here knows about
 * tetrahedra.
 *
 * Per-module imports throughout: the umbrella `@kitware/vtk.js` entry point
 * defeats tree-shaking and multiplies the bundle size.
 */
import '@kitware/vtk.js/Rendering/Profiles/Geometry';

import vtkActor from '@kitware/vtk.js/Rendering/Core/Actor';
import vtkMapper from '@kitware/vtk.js/Rendering/Core/Mapper';
import type vtkColorTransferFunction from '@kitware/vtk.js/Rendering/Core/ColorTransferFunction';
import type vtkPolyData from '@kitware/vtk.js/Common/DataModel/PolyData';
import vtkGenericRenderWindow from '@kitware/vtk.js/Rendering/Misc/GenericRenderWindow';
import vtkXMLPolyDataReader from '@kitware/vtk.js/IO/XML/XMLPolyDataReader';

import type { ArrayEntry, ScalarRange, Vector3 } from '../types';
import { DEFAULT_COLORMAP, makeColorTransferFunction } from './colormaps';
import {
    SOLID_COLOR,
    arrayOf,
    colorKey,
    computeRange,
    configureScalarMapper,
    listArrays,
} from './scalars';

export { SOLID_COLOR, colorKey };

/** How the surface is drawn. Matches vtk.js's `Property.Representation`. */
export const Representation = { POINTS: 0, WIREFRAME: 1, SURFACE: 2 } as const;
export type RepresentationMode = (typeof Representation)[keyof typeof Representation];

/** Default body colour when no array is selected. */
const BODY_COLOR: Vector3 = [0.85, 0.87, 0.9];

export interface LoadResult {
    numPoints: number;
    numCells: number;
    arrays: ArrayEntry[];
}

export class Renderer {
    private readonly window: ReturnType<typeof vtkGenericRenderWindow.newInstance>;
    private readonly renderer: ReturnType<
        ReturnType<typeof vtkGenericRenderWindow.newInstance>['getRenderer']
    >;
    private readonly renderWindow: ReturnType<
        ReturnType<typeof vtkGenericRenderWindow.newInstance>['getRenderWindow']
    >;

    private readonly bodyMapper: vtkMapper;
    private readonly bodyActor: vtkActor;
    /**
     * Edges are a separate actor, not `setEdgeVisibility` on the body.
     *
     * The mesh keeps its quads (`convertSurface` linearizes but does not
     * triangulate), and edge visibility draws the edges of the *triangulated*
     * primitives — which puts a fan diagonal across every quad.
     */
    private readonly edgeMapper: vtkMapper;
    private readonly edgeActor: vtkActor;

    private lut: vtkColorTransferFunction;
    private colormapName = DEFAULT_COLORMAP;
    private polydata: vtkPolyData | null = null;
    private currentEntry: ArrayEntry | null = null;
    private currentRange: ScalarRange | null = null;

    constructor(container: HTMLElement) {
        this.window = vtkGenericRenderWindow.newInstance({
            background: [0.09, 0.1, 0.13],
        });
        this.window.setContainer(container);
        this.renderer = this.window.getRenderer();
        this.renderWindow = this.window.getRenderWindow();

        this.lut = makeColorTransferFunction(this.colormapName, 0, 1);

        this.bodyMapper = vtkMapper.newInstance({
            useLookupTableScalarRange: true,
            scalarVisibility: false,
        });
        this.bodyMapper.setLookupTable(this.lut);

        this.bodyActor = vtkActor.newInstance({ mapper: this.bodyMapper });
        this.bodyActor.getProperty().setColor(...BODY_COLOR);
        this.renderer.addActor(this.bodyActor);

        this.edgeMapper = vtkMapper.newInstance({ scalarVisibility: false });
        // The edges lie exactly on the surface, so without an offset they
        // z-fight into a dashed mess. A *negative* offset moves toward the
        // viewer, so it belongs on the edges -- putting it on the body instead
        // pulls the surface in front of its own wireframe and hides it.
        this.edgeMapper.setResolveCoincidentTopologyToPolygonOffset();
        this.edgeMapper.setRelativeCoincidentTopologyLineOffsetParameters(-4, -4);
        this.edgeActor = vtkActor.newInstance({ mapper: this.edgeMapper });
        this.edgeActor.getProperty().setRepresentation(Representation.WIREFRAME);
        this.edgeActor.getProperty().setLineWidth(1);
        this.edgeActor.setVisibility(false);
        this.setEdgeColorFromBody();
        this.renderer.addActor(this.edgeActor);

        this.window.resize();
        window.addEventListener('resize', () => this.window.resize());
    }

    /** The render window vtk.js talks to WebGL through (for screenshots). */
    get view() {
        return this.window.getApiSpecificRenderWindow();
    }

    /** The event source widgets attach to. Lives on the generic window, not
     * on the API-specific one. */
    get interactor() {
        return this.window.getInteractor();
    }

    get scene() {
        return this.renderer;
    }

    get data(): vtkPolyData | null {
        return this.polydata;
    }

    /**
     * World-space extent of what is displayed.
     *
     * A boundary surface spans the same box as the solid it came from, so this
     * is a sound stand-in when a format's fast metadata path reports no
     * bounding box -- which VTU's does not.
     */
    bounds(): { min: Vector3; max: Vector3 } | null {
        if (!this.polydata) return null;
        const b = this.polydata.getBounds();
        if (!b || b.length < 6 || !Number.isFinite(b[0]) || b[0] > b[1]) return null;
        return { min: [b[0], b[2], b[4]], max: [b[1], b[3], b[5]] };
    }

    /** Replace the displayed mesh. */
    load(vtp: ArrayBuffer): LoadResult {
        const reader = vtkXMLPolyDataReader.newInstance();
        reader.parseAsArrayBuffer(vtp);
        const polydata = reader.getOutputData(0) as vtkPolyData;
        this.polydata = polydata;
        this.bodyMapper.setInputData(polydata);
        this.edgeMapper.setInputData(polydata);
        this.setColorBy(SOLID_COLOR);
        this.resetCamera();
        return {
            numPoints: polydata.getNumberOfPoints(),
            numCells: polydata.getNumberOfCells(),
            arrays: this.arrays(),
        };
    }

    arrays(): ArrayEntry[] {
        return listArrays(this.polydata);
    }

    /**
     * Colour by one entry of {@link arrays}, or `SOLID_COLOR` for a plain
     * surface. Returns the range that was applied, for the legend.
     */
    setColorBy(key: string): ScalarRange | null {
        const entry = this.arrays().find((a) => colorKey(a) === key) ?? null;
        this.currentEntry = entry;
        if (!entry) {
            this.currentRange = null;
            this.bodyMapper.setScalarVisibility(false);
            this.render();
            return null;
        }

        const array = arrayOf(this.polydata, entry);
        if (!array) {
            this.currentEntry = null;
            this.currentRange = null;
            this.bodyMapper.setScalarVisibility(false);
            this.render();
            return null;
        }

        // Range first, mapper second. `interpolateScalarsBeforeMapping`
        // caches a CPU-mapped colour array, so configuring the mapper against
        // a stale range renders the mesh in the wrong half of the colormap --
        // with no error and no obvious symptom beyond "the colours look odd".
        const range = computeRange(array, entry.component);
        this.currentRange = range;
        this.setLutRange(range.min, range.max);
        configureScalarMapper(this.bodyMapper, this.lut, entry);
        this.render();
        return range;
    }

    /** Narrow or widen the colour range without changing the array. */
    applyRange(min: number, max: number): void {
        this.setLutRange(min, max);
        // The cached colour array is keyed on the mapper, not the LUT, so it
        // has to be invalidated explicitly when only the range changes.
        this.bodyMapper.modified();
        this.render();
    }

    private setLutRange(min: number, max: number): void {
        // A degenerate range leaves the LUT with nothing to interpolate over.
        const hi = max > min ? max : min + 1e-12;
        this.lut.setMappingRange(min, hi);
        this.lut.updateRange();
    }

    /** The full data range of the current array, for "rescale to data". */
    dataRange(): ScalarRange | null {
        if (!this.currentEntry) return null;
        const array = arrayOf(this.polydata, this.currentEntry);
        return array ? computeRange(array, this.currentEntry.component) : null;
    }

    get range(): ScalarRange | null {
        return this.currentRange;
    }

    get colormap(): string {
        return this.colormapName;
    }

    setColormap(name: string): void {
        this.colormapName = name;
        const previous = this.lut.getMappingRange();
        this.lut = makeColorTransferFunction(name, previous[0], previous[1]);
        this.bodyMapper.setLookupTable(this.lut);
        if (this.currentEntry) {
            configureScalarMapper(this.bodyMapper, this.lut, this.currentEntry);
        }
        this.lut.setMappingRange(previous[0], previous[1]);
        this.lut.updateRange();
        this.render();
    }

    get lookupTable(): vtkColorTransferFunction {
        return this.lut;
    }

    setEdgeVisibility(visible: boolean): void {
        this.edgeActor.setVisibility(visible);
        this.render();
    }

    setRepresentation(mode: RepresentationMode): void {
        this.bodyActor.getProperty().setRepresentation(mode);
        // Edges over a wireframe or point cloud are noise, not information.
        this.edgeActor
            .getProperty()
            .setOpacity(mode === Representation.SURFACE ? 1 : 0);
        this.render();
    }

    setOpacity(opacity: number): void {
        this.bodyActor.getProperty().setOpacity(opacity);
        this.render();
    }

    setBackground(color: Vector3): void {
        this.renderer.setBackground(...color);
        this.setEdgeColorFromBody();
        this.render();
    }

    resetCamera(): void {
        // A three-quarter view, not vtk.js's default axis-aligned one: looking
        // straight down an axis renders a box as a flat rectangle, which reads
        // as a broken viewer rather than as a camera choice.
        const camera = this.renderer.getActiveCamera();
        camera.setPosition(1, -1, 0.6);
        camera.setFocalPoint(0, 0, 0);
        camera.setViewUp(0, 0, 1);
        this.renderer.resetCamera();
        this.renderer.resetCameraClippingRange();
        this.render();
    }

    render(): void {
        this.renderWindow.render();
    }

    resize(): void {
        this.window.resize();
    }

    private setEdgeColorFromBody(): void {
        // Half the body colour: bright enough to read against the surface,
        // dark enough not to compete with the data colours.
        // Dark enough not to compete with the data colours, bright enough to
        // read against them.
        const [r, g, b] = BODY_COLOR;
        this.edgeActor.getProperty().setColor(r * 0.22, g * 0.22, b * 0.22);
    }
}
