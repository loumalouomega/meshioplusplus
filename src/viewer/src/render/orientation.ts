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

/** Muted so the cube reads as a control, not as part of the mesh. */
const FACE_BASE = '#242833';

export interface OrientationHandle {
    /** Snap the camera along a world direction, as a face click would. */
    snapTo(normal: Vector3): void;
    setEnabled(enabled: boolean): void;
}

export function setupOrientationCube(
    renderer: Renderer,
    canvas: HTMLCanvasElement
): OrientationHandle {
    const scene = renderer.scene;
    const interactor = renderer.interactor;

    const cube = vtkAnnotatedCubeActor.newInstance();
    cube.setDefaultStyle({
        text: '',
        fontStyle: 'bold',
        fontFamily: 'Arial',
        fontColor: 'white',
        faceColor: FACE_BASE,
        // A dark groove between faces so adjacent faces stay distinguishable
        // when the cube is nearly edge-on.
        edgeThickness: 0.08,
        edgeColor: '#0b0d12',
        resolution: 400,
    });

    // Z-up: +Z is TOP, +Y is BACK. Axis colours follow the usual X=red,
    // Y=green, Z=blue convention so the cube and the arrows agree.
    cube.setXPlusFaceProperty({ text: 'RIGHT', faceColor: '#7a2222' });
    cube.setXMinusFaceProperty({ text: 'LEFT', faceColor: '#4a1414' });
    cube.setYPlusFaceProperty({ text: 'BACK', faceColor: '#1f6b2a' });
    cube.setYMinusFaceProperty({ text: 'FRONT', faceColor: '#134018' });
    cube.setZPlusFaceProperty({ text: 'TOP', faceColor: '#22447f' });
    cube.setZMinusFaceProperty({ text: 'BOTTOM', faceColor: '#122654' });

    const widget = vtkOrientationMarkerWidget.newInstance();
    widget.setActor(cube);
    widget.setInteractor(interactor);
    widget.setViewportCorner(vtkOrientationMarkerWidget.Corners.BOTTOM_LEFT);
    widget.setViewportSize(0.14);
    widget.setMinPixelSize(80);
    widget.setMaxPixelSize(160);
    widget.setEnabled(true);

    // Axes go into the *widget's* renderer, not the scene, so the arrows
    // rotate with the cube instead of sitting in the mesh.
    const axes = vtkAxesActor.newInstance();
    setAxesConfig(axes, {
        recenter: false,
        xLabel: 'X',
        yLabel: 'Y',
        zLabel: 'Z',
        tipLength: 0.25,
        tipRadius: 0.1,
        shaftRadius: 0.03,
    });
    // Anchored at the back-bottom-left corner and over-scaled so the tips
    // protrude past the opposing faces and stay readable head-on.
    placeAxes(axes, [-0.5, -0.5, -0.5], 1.65);
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

            picker.pick([displayX, displayY, 0], widget.getRenderer());
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
