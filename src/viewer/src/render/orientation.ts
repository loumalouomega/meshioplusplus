/**
 * The orientation cube: a labelled cube in the corner that shows which way is
 * up, and snaps the camera to an axis view when a face is clicked.
 *
 * vtk.js ships the widget but not the click-to-snap behaviour, so the hit-test
 * is done by hand. Three details in it are load-bearing and were established
 * the hard way in the sibling MDPA-Preview viewer:
 *
 * 1. The listener is on the **capture** phase and calls
 *    `stopImmediatePropagation()`, so clicking a cube face does not also start
 *    a camera rotate in the main scene.
 * 2. Display coordinates are **y-up**, browser coordinates are y-down, so the
 *    y must be flipped before picking.
 * 3. `widget.computeViewport()` returns `[left, bottom, right, top]` in 0..1,
 *    which is the only reliable way to know whether a click landed on the
 *    widget rather than the mesh.
 *
 * meshio++ is **Z-up** (see `Renderer.resetCamera`), unlike MDPA's Kratos Y-up
 * convention, so the face labels and the `viewUp` degeneracy guard differ.
 */
import vtkAnnotatedCubeActor from '@kitware/vtk.js/Rendering/Core/AnnotatedCubeActor';
import vtkAxesActor from '@kitware/vtk.js/Rendering/Core/AxesActor';
import vtkCellPicker from '@kitware/vtk.js/Rendering/Core/CellPicker';
import vtkOrientationMarkerWidget from '@kitware/vtk.js/Interaction/Widgets/OrientationMarkerWidget';

import type { Vector3 } from '../types';
import type { Renderer } from './renderer';
import { placeAxes, setAxesConfig } from '../vtk-shims';

/*
 * Look and colours follow the Kratos MDPA Preview extension's cube (itself
 * matched to CAD-Preview's), so the family reads alike: uniform light-blue
 * faces drawn unlit so the scene light cannot shade them unevenly, a darker
 * blue groove between faces, bold white labels sized to fit "BOTTOM", and the
 * conventional red/green/blue axis arrows emerging through the face centres.
 */
const FACE_COLOR = '#85b5da';
const EDGE_COLOR = '#5a87ae';
/** X red, Y green, Z blue, 0..255. */
const AXIS_COLORS: readonly [Vector3, Vector3, Vector3] = [
    [255, 54, 83],
    [138, 219, 0],
    [44, 143, 255],
];

export interface OrientationHandle {
    /** Snap the camera along a world direction, as a face click would. */
    snapTo(normal: Vector3): void;
    setEnabled(enabled: boolean): void;
}

export function setupOrientationCube(
    renderer: Renderer,
    canvas: HTMLCanvasElement,
    onSnap?: () => void
): OrientationHandle {
    const scene = renderer.scene;
    const interactor = renderer.interactor;

    const cube = vtkAnnotatedCubeActor.newInstance();
    cube.setDefaultStyle({
        text: '',
        fontStyle: 'bold',
        fontFamily: 'Arial',
        fontColor: 'white',
        // vtk.js's default scale overflows the face and crops longer words.
        fontSizeScale: (resolution: number) => resolution / 5,
        faceColor: FACE_COLOR,
        edgeThickness: 0.08,
        edgeColor: EDGE_COLOR,
        resolution: 400,
    });
    // Unlit: the faces show their true colour whatever the scene light does.
    cube.getProperty().setAmbient(1);
    cube.getProperty().setDiffuse(0);

    // Z-up: +Z is TOP, -Y is FRONT.
    cube.setXPlusFaceProperty({ text: 'RIGHT' });
    cube.setXMinusFaceProperty({ text: 'LEFT' });
    cube.setYPlusFaceProperty({ text: 'BACK' });
    cube.setYMinusFaceProperty({ text: 'FRONT' });
    cube.setZPlusFaceProperty({ text: 'TOP' });
    cube.setZMinusFaceProperty({ text: 'BOTTOM' });

    const widget = vtkOrientationMarkerWidget.newInstance();
    widget.setActor(cube);
    widget.setInteractor(interactor);
    widget.setViewportCorner(vtkOrientationMarkerWidget.Corners.BOTTOM_LEFT);
    widget.setViewportSize(0.15);
    widget.setMinPixelSize(80);
    widget.setMaxPixelSize(160);
    widget.setEnabled(true);

    // Axes go into the *widget's* renderer, not the scene, so the arrows
    // rotate with the cube. Anchored at its centre, each arrow leaves through
    // the middle of its face; the labelled faces already name the directions,
    // so the arrows carry no letters.
    const axes = vtkAxesActor.newInstance();
    setAxesConfig(axes, {
        recenter: false,
        xLabel: '',
        yLabel: '',
        zLabel: '',
        tipLength: 0.3,
        tipRadius: 0.12,
        shaftRadius: 0.035,
    });
    axes.setXAxisColor(AXIS_COLORS[0]);
    axes.setYAxisColor(AXIS_COLORS[1]);
    axes.setZAxisColor(AXIS_COLORS[2]);
    placeAxes(axes, [0, 0, 0], 1.15);
    widget.getRenderer().addActor(axes);

    const picker = vtkCellPicker.newInstance();

    function snapTo(normal: Vector3): void {
        const camera = scene.getActiveCamera();
        const focal = camera.getFocalPoint();
        const distance = camera.getDistance();
        camera.setPosition(
            focal[0] + normal[0] * distance,
            focal[1] + normal[1] * distance,
            focal[2] + normal[2] * distance
        );
        // Looking along ±Z makes the default Z viewUp parallel to the view
        // direction, which leaves the camera with no defined orientation.
        if (Math.abs(normal[2]) > 0.9) {
            camera.setViewUp(0, normal[2] > 0 ? 1 : -1, 0);
        } else {
            camera.setViewUp(0, 0, 1);
        }
        scene.resetCameraClippingRange();
        renderer.render();
        onSnap?.();
    }

    canvas.addEventListener(
        'pointerdown',
        (ev: PointerEvent) => {
            const rect = canvas.getBoundingClientRect();
            const displayX = ev.clientX - rect.left;
            // Display space is y-up; the browser gives y-down.
            const displayY = rect.height - (ev.clientY - rect.top);

            const vp = widget.computeViewport();
            const xNorm = displayX / rect.width;
            const yNorm = displayY / rect.height;
            const inWidget =
                xNorm >= vp[0] && xNorm <= vp[2] && yNorm >= vp[1] && yNorm <= vp[3];
            if (!inWidget) return;

            // Capture phase: stop vtk.js starting a rotate behind the widget.
            ev.stopImmediatePropagation();

            // The picker works in canvas pixels (CSS pixels x device pixel ratio).
            const sx = rect.width > 0 ? canvas.width / rect.width : 1;
            const sy = rect.height > 0 ? canvas.height / rect.height : 1;
            picker.pick([displayX * sx, displayY * sy, 0], widget.getRenderer());
            if (picker.getActors().length === 0) return;

            const n = picker.getMapperNormal();
            const length = Math.hypot(n[0], n[1], n[2]);
            if (length > 0.5) snapTo([n[0], n[1], n[2]]);
        },
        true
    );

    return {
        snapTo,
        setEnabled(enabled: boolean) {
            widget.setEnabled(enabled);
        },
    };
}
