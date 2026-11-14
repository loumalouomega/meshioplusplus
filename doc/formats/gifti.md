# GIFTI — neuroimaging surfaces (`.gii`)

The [GIFTI](https://www.nitrc.org/projects/gifti/) XML format for cortical and other surfaces and their per-vertex data, what FreeSurfer, Connectome Workbench, SimNIBS, Brainstorm and MNE exchange (roadmap 2.4.4). Both engines read and write it (`formats/gifti.cpp` through pugixml, `meshioplusplus.gifti`), and they write the same bytes.

| | |
|---|---|
| **Format name** | `gifti` |
| **Extensions** | `.gii` (`.surf.gii`, `.func.gii`, `.label.gii`, `.shape.gii` all end in it) |
| **Read / Write** | ✓ / ✓ (one triangle surface) |
| **Extra dependencies** | — (zlib for `GZipBase64Binary` arrays) |

## Reading

```python
import meshioplusplus

mesh = meshioplusplus.read("lh.pial.surf.gii")
```

Every `DataArray` is decoded: `ASCII`, `Base64Binary` and `GZipBase64Binary` (base64, then a zlib stream), little- or big-endian, row- or column-major, any NIfTI integer or float type. External-file arrays are refused by name.

- `NIFTI_INTENT_POINTSET` is the points; its `CoordinateSystemTransformMatrix` is `field_data["gifti:transform"]` (4x4) and its space names `gifti:data_space` and `gifti:transformed_space` (bytes).
- `NIFTI_INTENT_TRIANGLE` is the triangle block.
- Every other array with one row per vertex is point data named by its `Name` metadata, else `gifti:<intent>_<index>` (`gifti:label` for an unnamed `NIFTI_INTENT_LABEL`). An array of another length is skipped with a warning.
- With a `NIFTI_INTENT_LABEL` array, each `LabelTable` entry is a point region over the vertices carrying its key, tagged with the key.

## Writing

```python
meshioplusplus.write("out.surf.gii", mesh)
```

One triangle block (or none): the points as a float32 `POINTSET` with the transform (identity when the mesh has none), the triangles as an int32 `TRIANGLE` array, then each point data array as float32 (`SHAPE`, or `NONE` with several components) or int32 (`NONE`; `LABEL` for the one named `gifti:label`), every array little-endian `Base64Binary`, row-major and carrying its name. Point regions with a non-negative tag become the `LabelTable`. Coordinates and floating-point data are float32, as the format stores them.

## Verification

nibabel 5.5 reads what meshio++ writes (points, triangles, named arrays, the label table) and meshio++ reads what nibabel writes in every encoding and byte order (`tests/python/test_gifti.py`).
