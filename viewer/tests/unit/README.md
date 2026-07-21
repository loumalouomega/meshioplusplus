# Unit tests

Pure-logic tests for modules with no DOM and no WebGL — currently the colour
maps, whose whole point is that the legend and the 3D render are driven by the
same stops.

They run on `node --test`, which cannot import TypeScript directly until Node
22.6; this repo's CI is on Node 20. So `npm run test:unit` bundles each test
(and the `.ts` it imports) with esbuild — already present as a Vite dependency,
no new package — and runs the plain JS output.
