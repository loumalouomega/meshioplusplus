# COMSOL binary mesh (`.mphbin`)

The binary twin of the [COMSOL text mesh](./mphtxt.md): the same objects, element types, node order, entity indices and Selections, serialised in binary.

| | |
|---|---|
| **Format name** | `mphbin` |
| **Extensions** | `.mphbin` |
| **Read / Write** | ✓ / ✓ |
| **Extra dependencies** | — |

## Reading & writing

```python
import meshioplusplus

mesh = meshioplusplus.read("model.mphbin")
meshioplusplus.mphbin.write("out.mphbin", mesh)
```

## Encoding

- Integers are little-endian int32 and doubles little-endian float64.
- A string is an int32 length followed by one int32 code point per character.
- A file opens with the int32 values `0 1` (the version), then the tag count, which is how `sniff_format` recognises it.

Everything else — objects, Mesh versions, Selections, node order, what the writer derives — is described on the [mphtxt page](./mphtxt.md).

The binary layout follows COMSOL's [binary file format](https://doc.comsol.com/6.0/doc/com.comsol.help.comsol/comsol_api_fileformats.50.22.html) and [AWS Palace](https://github.com/awslabs/palace)'s reader. No public COMSOL-written `.mphbin` file was found to test against; the fixture `tests/python/meshes/comsol/two_domains.mphbin` is written by `tools/gen_comsol_fixtures.py`, which encodes the same values as `two_domains.mphtxt` independently of meshio++.

## Quirks & limitations

- **No provenance.** The format has no comment, so none is written. The writer still consults the provenance machinery, so `Mode.REQUIRED` raises there, as for the other formats without a slot.
- **32-bit values.** Values beyond 32 bits cannot be written.
