// A Kuhn-tetrahedralised unit cube as a JS mesh object, shared by the WebAssembly
// benchmark (bench.mjs) and the memory test (memory.mjs).

// The Kuhn tetrahedralisation of the unit cube: 6 * n^3 tetrahedra.
export function kuhnCube(n) {
    const np = n + 1;
    const points = new Float64Array(np * np * np * 3);
    for (let k = 0, p = 0; k < np; k++)
        for (let j = 0; j < np; j++)
            for (let i = 0; i < np; i++, p += 3) {
                points[p] = i / n;
                points[p + 1] = j / n;
                points[p + 2] = k / n;
            }
    const id = (i, j, k) => (k * np + j) * np + i;
    const paths = [
        [[1, 0, 0], [1, 1, 0]], [[1, 0, 0], [1, 0, 1]], [[0, 1, 0], [1, 1, 0]],
        [[0, 1, 0], [0, 1, 1]], [[0, 0, 1], [1, 0, 1]], [[0, 0, 1], [0, 1, 1]],
    ];
    const conn = new Int32Array(6 * n * n * n * 4);
    let c = 0;
    for (let k = 0; k < n; k++)
        for (let j = 0; j < n; j++)
            for (let i = 0; i < n; i++)
                for (const [a, b] of paths) {
                    const v = [[0, 0, 0], a, b, [1, 1, 1]].map(([x, y, z]) => id(i + x, j + y, k + z));
                    // Orient every tetrahedron positively, so the file is a valid mesh.
                    const p = v.map((q) => points.subarray(q * 3, q * 3 + 3));
                    const d = [1, 2, 3].map((m) => [0, 1, 2].map((x) => p[m][x] - p[0][x]));
                    const det =
                        d[0][0] * (d[1][1] * d[2][2] - d[1][2] * d[2][1]) -
                        d[0][1] * (d[1][0] * d[2][2] - d[1][2] * d[2][0]) +
                        d[0][2] * (d[1][0] * d[2][1] - d[1][1] * d[2][0]);
                    if (det < 0) [v[2], v[3]] = [v[3], v[2]];
                    conn.set(v, c);
                    c += 4;
                }
    const nTet = c / 4;
    const u = new Float64Array(np * np * np);
    for (let q = 0; q < u.length; q++) u[q] = points[q * 3] ** 2 + points[q * 3 + 1] * points[q * 3 + 2];
    const material = new Float64Array(nTet).map((_, q) => q % 7);
    return {
        points,
        dim: 3,
        cells: [{ type: 'tetra', data: conn, nodesPerCell: 4 }],
        point_data: { u },
        cell_data: { material: [material] },
        field_data: {},
    };
}
