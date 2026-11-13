# FreeSurfer — binary surfaces, curv and annot files

FreeSurfer's binary triangle surfaces (`lh.pial`, `rh.white`, `lh.inflated`, `lh.sphere`, …) with their `curv` overlays and `annot` parcellations (roadmap 2.4.4). Both engines read and write the surface (`formats/freesurfer.cpp`, `meshioplusplus.freesurfer`).

| | |
|---|---|
| **Format name** | `freesurfer` |
| **Extensions** | `.pial`, `.white`, `.inflated`, `.orig`, `.smoothwm`, `.sphere` (FreeSurfer names surfaces by their last suffix) |
| **Read / Write** | ✓ (surface, curv, annot) / ✓ (surface) |
| **Extra dependencies** | — |

## Reading

```python
import meshioplusplus

mesh = meshioplusplus.read("surf/lh.pial")
mesh = meshioplusplus.freesurfer.read(
    "surf/lh.white",
    curv=["surf/lh.thickness", "surf/lh.curv"],
    annot="label/lh.aparc.annot",
)
```

- The surface: the magic `0xFFFFFE`, a creator line ended by two newlines (`field_data["freesurfer:created_by"]`), big-endian vertex and face counts, float32 coordinates and int32 faces, as one triangle block. Whatever follows the faces (the volume-geometry tags FreeSurfer appends) is kept verbatim in `field_data["freesurfer:tail"]`. Quad surfaces are refused by name.
- `curv=`: each new-format curv file (`0xFFFFFF`, one value per vertex) is float32 point data named after the file without its `lh.`/`rh.` prefix (`lh.thickness` is `thickness`). The old format is refused by name.
- `annot=`: each vertex's colour is matched to the colour table (either layout: the old one, or version 2), giving the Int32 point data `freesurfer:annot` (the structure index, -1 for a colour no entry has) and one point region per entry, tagged with its index.

## Writing

```python
meshioplusplus.write("lh.out.pial", mesh)
```

One triangle block, with the creator line and the kept tail when the mesh has them. A surface read from a FreeSurfer file writes back byte-identical. Overlays are not written.

## Verification

nibabel 5.5 writes the surface, curv and annot files the tests read, and reads the surface meshio++ writes, byte for byte (`tests/python/test_freesurfer.py`). A real `fsaverage` subject is under the FreeSurfer licence and is not copied here.
