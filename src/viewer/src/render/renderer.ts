/**
 * The vtk.js scene.
 *
 * vtk.js's only mesh data model is `vtkPolyData` — it ships no unstructured
 * grid and no XML unstructured-grid reader — which is why a volume mesh is
 * shown by its boundary surface. The meshio++ WASM side does that extraction
 * and hands the surface over as the arrays a `vtkPolyData` is built from
 * ({@link Renderer.loadSurface}); the offline embedded page still carries a
 * VTP ({@link Renderer.load}). Nothing here knows about tetrahedra.
 *
 * Per-module imports throughout: the umbrella `@kitware/vtk.js` entry point
 * defeats tree-shaking and multiplies the bundle size.
 */
import '@kitware/vtk.js/Rendering/Profiles/Geometry';

import vtkActor from '@kitware/vtk.js/Rendering/Core/Actor';
import vtkMapper from '@kitware/vtk.js/Rendering/Core/Mapper';
import vtkColorTransferFunction from '@kitware/vtk.js/Rendering/Core/ColorTransferFunction';
import vtkCellArray from '@kitware/vtk.js/Common/Core/CellArray';
import vtkDataArray from '@kitware/vtk.js/Common/Core/DataArray';
import vtkPlane from '@kitware/vtk.js/Common/DataModel/Plane';
import vtkPolyData from '@kitware/vtk.js/Common/DataModel/PolyData';
import { Scale } from '@kitware/vtk.js/Rendering/Core/ColorTransferFunction/Constants';
import vtkGenericRenderWindow from '@kitware/vtk.js/Rendering/Misc/GenericRenderWindow';
import vtkXMLPolyDataReader from '@kitware/vtk.js/IO/XML/XMLPolyDataReader';

import type { SurfaceBuffers } from '@meshioplusplus/wasm';

import type { ArrayEntry, ScalarRange, Vector3 } from '../types';
import type { CellTopology, PointCells } from '../selection/topology';
import {
    buildCellTopology,
    buildPointCells,
    cellCentroids,
    cellNormals,
} from '../selection/topology';
import { canLogScale } from './scale';
import { subsetCells, takeTuples } from './subset';
import { type CameraState, type ViewName, viewPose } from './camera';
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

/** The per-cell label array behind colour-by-group (never offered in the colour-by menu). */
const LABEL_ARRAY = 'viewer:label';

/** How the surface is drawn. Matches vtk.js's `Property.Representation`. */
export const Representation = { POINTS: 0, WIREFRAME: 1, SURFACE: 2 } as const;
export type RepresentationMode = (typeof Representation)[keyof typeof Representation];

/** Default body colour when no array is selected (until a theme sets one). */
const BODY_COLOR: Vector3 = [0.78, 0.82, 0.88];

/** A live clipping plane: an axis, a position along it (0..1 over the bounds) and a side. */
export interface ClipSpec {
    axis: 0 | 1 | 2;
    t: number;
    flip: boolean;
}

export interface LightingSpec {
    ambient: number;
    diffuse: number;
    specular: number;
    specularPower: number;
}

/** What screen-space selection needs from the camera, in CSS pixels. */
export interface CameraProjection {
    /** Row-major world -> clip matrix. */
    matrix: number[];
    width: number;
    height: number;
    eye: Vector3;
    /** Unit view direction (the way the camera looks). */
    direction: Vector3;
    parallel: boolean;
}

/** Per-surface geometry the selection tools use, built on first need. */
export interface SurfaceGeometry {
    topo: CellTopology;
    centroids: Float64Array;
    normals: Float32Array;
    pointCells: PointCells;
}

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
    /** The full surface. Ranges, picking ids and exports always refer to this. */
    private polydata: vtkPolyData | null = null;
    /** What the mappers draw: `polydata`, or a subset of it while cells are hidden. */
    private shown: vtkPolyData | null = null;
    private cellMap: Uint32Array | null = null;
    private visibleMask: Uint8Array | null = null;
    private scaleLog = false;
    private scaleBands = 0;
    /** Original point coordinates while a deformation is applied. */
    private restPoints: Float32Array | null = null;
    private readonly clipPlane = vtkPlane.newInstance();
    private clip: ClipSpec | null = null;
    private readonly clippables = new Set<vtkMapper>();
    private readonly geometryListeners = new Set<() => void>();
    private geometry: SurfaceGeometry | null = null;
    private currentEntry: ArrayEntry | null = null;
    private currentRange: ScalarRange | null = null;
    private bodyColor: Vector3 = BODY_COLOR;
    private background: Vector3 = [0.118, 0.118, 0.118];

    constructor(container: HTMLElement) {
        this.window = vtkGenericRenderWindow.newInstance({
            background: [0.118, 0.118, 0.118],
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
        this.bodyActor.getProperty().setColor(...this.bodyColor);
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

        this.makeClippable(this.bodyMapper);
        this.makeClippable(this.edgeMapper);

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

    /** Replace the displayed mesh with a VTP file's (the offline embedded page). */
    load(vtp: ArrayBuffer): LoadResult {
        const reader = vtkXMLPolyDataReader.newInstance();
        reader.parseAsArrayBuffer(vtp);
        return this.show(reader.getOutputData(0) as vtkPolyData);
    }

    /**
     * Replace the displayed mesh with a surface the worker built. The arrays
     * are laid out as vtk.js's own VTP reader lays them out, so this and
     * {@link load} produce the same `vtkPolyData` for the same surface.
     */
    loadSurface(surface: SurfaceBuffers): LoadResult {
        const polydata = vtkPolyData.newInstance();
        polydata.getPoints().setData(surface.points, 3);
        if (surface.verts.length) polydata.setVerts(vtkCellArray.newInstance({ values: surface.verts }));
        if (surface.lines.length) polydata.setLines(vtkCellArray.newInstance({ values: surface.lines }));
        if (surface.polys.length) polydata.setPolys(vtkCellArray.newInstance({ values: surface.polys }));
        for (const a of surface.pointData) {
            polydata.getPointData().addArray(
                vtkDataArray.newInstance({ name: a.name, values: a.values, numberOfComponents: a.components })
            );
        }
        for (const a of surface.cellData) {
            polydata.getCellData().addArray(
                vtkDataArray.newInstance({ name: a.name, values: a.values, numberOfComponents: a.components })
            );
        }
        return this.show(polydata);
    }

    private show(polydata: vtkPolyData): LoadResult {
        this.polydata = polydata;
        this.shown = polydata;
        this.cellMap = null;
        this.visibleMask = null;
        this.geometry = null;
        this.restPoints = null;
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

    /** The raw VTK legacy cell arrays of the full surface. */
    private rawCells(): { verts: Uint32Array; lines: Uint32Array; polys: Uint32Array } | null {
        const p = this.polydata;
        if (!p) return null;
        const get = (cells: { getData(): ArrayLike<number> } | null | undefined) =>
            Uint32Array.from(cells ? cells.getData() : []);
        return { verts: get(p.getVerts()), lines: get(p.getLines()), polys: get(p.getPolys()) };
    }

    /** Topology, centroids and normals of the full surface (cached per surface). */
    surfaceGeometry(): SurfaceGeometry | null {
        if (this.geometry) return this.geometry;
        const p = this.polydata;
        const raw = this.rawCells();
        if (!p || !raw) return null;
        const coords = p.getPoints().getData();
        const topo = buildCellTopology(raw.verts, raw.lines, raw.polys);
        this.geometry = {
            topo,
            centroids: cellCentroids(coords, topo),
            normals: cellNormals(coords, topo),
            pointCells: buildPointCells(topo, p.getNumberOfPoints()),
        };
        return this.geometry;
    }

    /**
     * Show only the cells in `mask` (1 = visible); `null` shows everything.
     * Points are shared with the full surface; per-cell arrays are restricted
     * to the subset so colouring by them still works.
     */
    setVisibleMask(mask: Uint8Array | null): void {
        this.visibleMask = mask;
        const full = this.polydata;
        const raw = this.rawCells();
        if (!full || !raw) return;
        if (!mask) {
            this.shown = full;
            this.cellMap = null;
        } else {
            const subset = subsetCells(raw.verts, raw.lines, raw.polys, mask);
            const sub = vtkPolyData.newInstance();
            sub.setPoints(full.getPoints());
            if (subset.verts.length) sub.setVerts(vtkCellArray.newInstance({ values: subset.verts }));
            if (subset.lines.length) sub.setLines(vtkCellArray.newInstance({ values: subset.lines }));
            if (subset.polys.length) sub.setPolys(vtkCellArray.newInstance({ values: subset.polys }));
            const pd = full.getPointData();
            for (let i = 0; i < pd.getNumberOfArrays(); i++) {
                const a = pd.getArrayByIndex(i);
                if (a) sub.getPointData().addArray(a);
            }
            const cd = full.getCellData();
            for (let i = 0; i < cd.getNumberOfArrays(); i++) {
                const a = cd.getArrayByIndex(i);
                if (!a) continue;
                const comps = a.getNumberOfComponents();
                sub.getCellData().addArray(
                    vtkDataArray.newInstance({
                        name: a.getName(),
                        numberOfComponents: comps,
                        values: takeTuples(a.getData() as Float32Array, comps, subset.cellMap),
                    })
                );
            }
            this.shown = sub;
            this.cellMap = subset.cellMap;
        }
        this.bodyMapper.setInputData(this.shown);
        this.edgeMapper.setInputData(this.shown);
        this.bodyMapper.modified();
        this.render();
    }

    /** Fit the camera to the cells in `mask` (1 = include), keeping the view direction. */
    frameCells(mask: Uint8Array): void {
        const g = this.surfaceGeometry();
        const p = this.polydata;
        if (!g || !p) return;
        const coords = p.getPoints().getData();
        const b: [number, number, number, number, number, number] = [
            Infinity,
            -Infinity,
            Infinity,
            -Infinity,
            Infinity,
            -Infinity,
        ];
        for (let c = 0; c < g.topo.numCells; c++) {
            if (!mask[c]) continue;
            for (let i = g.topo.offsets[c]; i < g.topo.offsets[c + 1]; i++) {
                const q = g.topo.points[i] * 3;
                for (let k = 0; k < 3; k++) {
                    b[k * 2] = Math.min(b[k * 2], coords[q + k]);
                    b[k * 2 + 1] = Math.max(b[k * 2 + 1], coords[q + k]);
                }
            }
        }
        if (!Number.isFinite(b[0])) return;
        this.renderer.resetCamera(b);
        this.renderer.resetCameraClippingRange();
        this.render();
    }

    /** The full-surface cell id behind a cell id of what is drawn. */
    toOriginalCell(id: number): number {
        return this.cellMap ? (this.cellMap[id] ?? -1) : id;
    }

    /** The camera as selection needs it, or null before anything is drawn. */
    cameraProjection(): CameraProjection | null {
        const canvas = this.view.getCanvas();
        if (!canvas) return null;
        const width = canvas.clientWidth;
        const height = canvas.clientHeight;
        if (width <= 0 || height <= 0) return null;
        const camera = this.renderer.getActiveCamera();
        const m = camera.getCompositeProjectionMatrix(width / height, -1, 1);
        const pos = camera.getPosition();
        const dir = camera.getDirectionOfProjection();
        return {
            matrix: Array.from(m),
            width,
            height,
            eye: [pos[0], pos[1], pos[2]],
            direction: [dir[0], dir[1], dir[2]],
            parallel: camera.getParallelProjection(),
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
        // A previous categorical colouring swapped the lookup table; the normal
        // colour maps always go through the shared one.
        this.bodyMapper.setLookupTable(this.lut);
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

    /**
     * Colour each cell by an integer label (group index), one distinct colour
     * per label; -1 draws the neutral grey. Replaces any colouring by an array.
     */
    setLabelColoring(labels: Int32Array, colors: Vector3[]): void {
        const full = this.polydata;
        if (!full) return;
        const values = new Float32Array(labels.length);
        for (let i = 0; i < labels.length; i++) values[i] = labels[i] < 0 ? Number.NaN : labels[i];
        full.getCellData().removeArray(LABEL_ARRAY);
        full.getCellData().addArray(
            vtkDataArray.newInstance({ name: LABEL_ARRAY, values, numberOfComponents: 1 })
        );
        // Rebuild what is drawn so the subset (if any) carries the new array.
        this.setVisibleMask(this.visibleMask);

        // A step function over integer labels: each label owns [i-0.5, i+0.5].
        const lut = vtkColorTransferFunction.newInstance();
        colors.forEach(([r, g, b], i) => {
            lut.addRGBPoint(i - 0.4999, r, g, b);
            lut.addRGBPoint(i + 0.4999, r, g, b);
        });
        lut.setNanColor(0.55, 0.55, 0.55, 1);
        lut.setMappingRange(-0.5, Math.max(0.5, colors.length - 0.5));
        lut.updateRange();
        this.bodyMapper.setLookupTable(lut);
        configureScalarMapper(this.bodyMapper, lut, {
            location: 'cell',
            name: LABEL_ARRAY,
            component: -1,
            label: '',
        });
        this.currentEntry = null;
        this.currentRange = null;
        this.bodyMapper.modified();
        this.render();
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
        this.applyScalarOptions();
    }

    /** The data range of any colourable array, without colouring by it. */
    rangeOf(entry: ArrayEntry): ScalarRange | null {
        const array = arrayOf(this.polydata, entry);
        return array ? computeRange(array, entry.component) : null;
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
        this.applyScalarOptions();
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
        this.background = color;
        this.renderer.setBackground(...color);
        this.setEdgeColorFromBody();
        this.render();
    }

    /** Scene colours for a colour scheme: background plus the plain-surface colour. */
    setSceneTheme(background: Vector3, body: Vector3): void {
        this.bodyColor = body;
        this.bodyActor.getProperty().setColor(...body);
        this.setBackground(background);
    }

    get backgroundColor(): Vector3 {
        return this.background;
    }

    resetCamera(view: ViewName = 'iso'): void {
        // The default is a three-quarter view, not vtk.js's axis-aligned one:
        // looking straight down an axis renders a box as a flat rectangle,
        // which reads as a broken viewer rather than as a camera choice.
        const pose = viewPose(view);
        const camera = this.renderer.getActiveCamera();
        camera.setPosition(...pose.direction);
        camera.setFocalPoint(0, 0, 0);
        camera.setViewUp(...pose.up);
        this.renderer.resetCamera();
        this.renderer.resetCameraClippingRange();
        this.render();
    }

    get parallel(): boolean {
        return this.renderer.getActiveCamera().getParallelProjection();
    }

    setParallel(on: boolean): void {
        const camera = this.renderer.getActiveCamera();
        if (camera.getParallelProjection() === on) return;
        if (on) {
            // Match the apparent size of what the perspective camera frames.
            const p = camera.getPosition();
            const f = camera.getFocalPoint();
            const distance = Math.hypot(p[0] - f[0], p[1] - f[1], p[2] - f[2]);
            const angle = (camera.getViewAngle() * Math.PI) / 180;
            camera.setParallelScale(distance * Math.tan(angle / 2));
        }
        camera.setParallelProjection(on);
        this.renderer.resetCameraClippingRange();
        this.render();
    }

    getCameraState(): CameraState {
        const c = this.renderer.getActiveCamera();
        const v = (a: number[]): [number, number, number] => [a[0], a[1], a[2]];
        return {
            position: v(c.getPosition()),
            focalPoint: v(c.getFocalPoint()),
            viewUp: v(c.getViewUp()),
            parallel: c.getParallelProjection(),
            parallelScale: c.getParallelScale(),
            viewAngle: c.getViewAngle(),
        };
    }

    setCameraState(state: CameraState): void {
        const c = this.renderer.getActiveCamera();
        c.setPosition(...state.position);
        c.setFocalPoint(...state.focalPoint);
        c.setViewUp(...state.viewUp);
        c.setParallelProjection(state.parallel);
        if (state.parallelScale) c.setParallelScale(state.parallelScale);
        if (state.viewAngle) c.setViewAngle(state.viewAngle);
        this.renderer.resetCameraClippingRange();
        this.render();
    }

    // --- colour scale -------------------------------------------------------

    /** Log scaling and banding of the colour map (vtk.js's own LOG10 and discretisation). */
    setScalarOptions(options: { log: boolean; bands: number }): void {
        this.scaleLog = options.log;
        this.scaleBands = options.bands;
        this.applyScalarOptions();
        this.bodyMapper.modified();
        this.render();
    }

    /** Whether the log scale is in effect: it needs a strictly positive range. */
    get logActive(): boolean {
        if (!this.scaleLog) return false;
        const [lo, hi] = this.lut.getMappingRange();
        return canLogScale(lo, hi);
    }

    private applyScalarOptions(): void {
        this.lut.setScale(this.logActive ? Scale.LOG10 : Scale.LINEAR);
        this.lut.setDiscretize(this.scaleBands > 0);
        if (this.scaleBands > 0) this.lut.setNumberOfValues(this.scaleBands);
    }

    // --- lighting -----------------------------------------------------------

    setLighting(spec: LightingSpec): void {
        const p = this.bodyActor.getProperty();
        p.setAmbient(spec.ambient);
        p.setDiffuse(spec.diffuse);
        p.setSpecular(spec.specular);
        p.setSpecularPower(spec.specularPower);
        this.render();
    }

    setBackfaceCulling(on: boolean): void {
        this.bodyActor.getProperty().setBackfaceCulling(on);
        this.render();
    }

    // --- deformation ----------------------------------------------------------

    /** Names of the point arrays that can displace the mesh (2 or 3 components). */
    vectorArrays(): string[] {
        const pd = this.polydata?.getPointData();
        if (!pd) return [];
        const out: string[] = [];
        for (let i = 0; i < pd.getNumberOfArrays(); i++) {
            const a = pd.getArrayByIndex(i);
            if (a && [2, 3].includes(a.getNumberOfComponents()) && !a.getName().startsWith('viewer:')) {
                out.push(a.getName());
            }
        }
        return out;
    }

    /**
     * Displace the points by a vector point array, scaled so that
     * `percent = 100` moves the largest displacement by 10% of the model's
     * diagonal (so the result is visible whatever the field's units are).
     * `name = null` puts the points back. Returns the largest displacement drawn.
     */
    setDeform(name: string | null, percent: number): number {
        const p = this.polydata;
        if (!p) return 0;
        const points = p.getPoints();
        if (!this.restPoints) this.restPoints = Float32Array.from(points.getData());
        const rest = this.restPoints;

        const array = name ? p.getPointData().getArrayByName(name) : null;
        let drawn = 0;
        if (!array) {
            points.setData(rest.slice(), 3);
        } else {
            const comps = array.getNumberOfComponents();
            const vec = array.getData();
            const n = Math.floor(rest.length / 3);
            let maxMag = 0;
            for (let i = 0; i < n; i++) {
                let s = 0;
                for (let k = 0; k < comps; k++) s += vec[i * comps + k] ** 2;
                maxMag = Math.max(maxMag, Math.sqrt(s));
            }
            const b = p.getBounds();
            const diag = Math.hypot(b[1] - b[0], b[3] - b[2], b[5] - b[4]);
            const scale = maxMag > 0 ? ((percent / 100) * 0.1 * diag) / maxMag : 0;
            const next = new Float32Array(rest.length);
            for (let i = 0; i < n; i++) {
                for (let k = 0; k < 3; k++) {
                    next[i * 3 + k] = rest[i * 3 + k] + (k < comps ? scale * vec[i * comps + k] : 0);
                }
            }
            points.setData(next, 3);
            drawn = maxMag * scale;
        }
        points.modified();
        p.modified();
        this.shown?.modified();
        this.geometry = null;
        this.bodyMapper.modified();
        this.edgeMapper.modified();
        this.renderer.resetCameraClippingRange();
        this.render();
        for (const listener of [...this.geometryListeners]) listener();
        return drawn;
    }

    /** Called after the point positions changed (a deformation was applied or removed). */
    onGeometryChanged(listener: () => void): () => void {
        this.geometryListeners.add(listener);
        return () => this.geometryListeners.delete(listener);
    }

    // --- clipping -----------------------------------------------------------------

    /** Make a mapper honour the live clip plane (now and whenever it changes). */
    makeClippable(mapper: vtkMapper): void {
        this.clippables.add(mapper);
        mapper.removeAllClippingPlanes();
        if (this.clip) mapper.addClippingPlane(this.clipPlane);
    }

    forgetClippable(mapper: vtkMapper): void {
        this.clippables.delete(mapper);
    }

    /** Clip the surface by a plane perpendicular to an axis; `null` turns clipping off. */
    setClip(spec: ClipSpec | null): void {
        const wasOn = this.clip !== null;
        this.clip = spec;
        if (spec) {
            const b = this.polydata?.getBounds() ?? [0, 1, 0, 1, 0, 1];
            const lo = b[spec.axis * 2];
            const hi = b[spec.axis * 2 + 1];
            const origin: [number, number, number] = [0, 0, 0];
            origin[spec.axis] = lo + Math.min(1, Math.max(0, spec.t)) * (hi - lo);
            const normal: [number, number, number] = [0, 0, 0];
            // vtk keeps what the normal points away from; flip swaps the side.
            normal[spec.axis] = spec.flip ? 1 : -1;
            this.clipPlane.setOrigin(...origin);
            this.clipPlane.setNormal(...normal);
        }
        if (wasOn !== (spec !== null)) {
            for (const mapper of this.clippables) {
                mapper.removeAllClippingPlanes();
                if (spec) mapper.addClippingPlane(this.clipPlane);
            }
        }
        for (const mapper of this.clippables) mapper.modified();
        this.render();
    }

    /** Whether a world position lies on the clipped-away side of the live plane. */
    isClipped(position: [number, number, number]): boolean {
        if (!this.clip) return false;
        const n = this.clipPlane.getNormal();
        const o = this.clipPlane.getOrigin();
        const d = n[0] * (position[0] - o[0]) + n[1] * (position[1] - o[1]) + n[2] * (position[2] - o[2]);
        return d > 0;
    }

    render(): void {
        this.renderWindow.render();
    }

    resize(): void {
        this.window.resize();
    }

    private setEdgeColorFromBody(): void {
        // A dark line over the surface colour reads against both a light and a
        // dark scene, and stays quieter than the data colours.
        const [r, g, b] = this.bodyColor;
        this.edgeActor.getProperty().setColor(r * 0.25, g * 0.25, b * 0.25);
    }
}
