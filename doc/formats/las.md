# LAS — ASPRS point clouds (`.las`)

The ASPRS LAS format airborne and mobile laser scanners write, versions 1.2 to 1.4 and point formats 0 to 10 (roadmap 2.4.6.1). Both engines read and write it natively (`formats/las.cpp`, `meshioplusplus.las`), and they write the same bytes. LAZ (compressed LAS) needs LASzip and is refused by name.

| | |
|---|---|
| **Format name** | `las` |
| **Extensions** | `.las` |
| **Read / Write** | ✓ / ✓ (LAS 1.4, point formats 6 to 8) |
| **Extra dependencies** | — |

## Reading

```python
import meshioplusplus

cloud = meshioplusplus.read("scan.las")
thin = meshioplusplus.las.read("scan.las", every=10, arrays=["classification", "intensity"])
```

- The points are `X Y Z` scaled by the header's scale and offset, in one `vertex` block (as [`pcd`](./pcd.md) does).
- Every field of the point format is point data under laspy's names and in its LAS type: `intensity`, `return_number`, `number_of_returns`, `scan_direction_flag`, `edge_of_flight_line`, `classification`, `synthetic`, `key_point`, `withheld`, `overlap` and `scanner_channel` (formats 6 to 10), `scan_angle_rank` (formats 0 to 5) or `scan_angle` (6 to 10, in 0.006-degree steps), `user_data`, `point_source_id`, `gps_time`, `red`, `green`, `blue` and `nir`. Wave-packet fields and any bytes past the standard record length are skipped.
- The variable-length records are kept verbatim in `field_data["las:vlrs"]` (bytes), so the CRS (GeoKeys or WKT) travels with the cloud; extended VLRs are not kept, with a warning. `las:scale`, `las:offset`, `las:version` and `las:point_format` are field data too.
- `every=n` keeps every n-th point (`LasReadOptions::mEvery`); `arrays=` and `points_only=` are the usual selective-read options.

## Writing

```python
meshioplusplus.write("out.las", cloud)
```

LAS 1.4 with point format 6, or 7 when the mesh has `red`, `green` and `blue`, or 8 with `nir` too. Every point of the mesh is written (its cells are not). The scale and offset come from `las:scale` and `las:offset`, else 0.001 and the floored minimum; a coordinate that does not fit the int32 grid is refused by name. The kept VLRs are written back, and the fields the mesh carries fill the record (`return_number` and `number_of_returns` default to 1, `scan_angle_rank` is converted to `scan_angle`). The creation date and GUID are zero, so the bytes depend on the cloud alone.

## Verification

laspy 2.7 writes the files the tests read (point formats 0, 1, 2, 3, 6, 7 and 8; versions 1.2 to 1.4), and reads what meshio++ writes (`tests/python/test_las.py`).
