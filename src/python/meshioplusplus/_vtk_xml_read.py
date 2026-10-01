"""Shared, private VTK XML framing for the PolyData and structured readers."""

from ._exceptions import ReadError
from .vtu._vtu import VtuReader, _compressor_for, _load_root


class ArrayReader:
    """Reuse VTU's decoder rather than maintaining another binary framing path."""

    _ordered = VtuReader._ordered
    read_appended = VtuReader.read_appended
    read_uncompressed_binary = VtuReader.read_uncompressed_binary
    read_compressed_binary = VtuReader.read_compressed_binary
    read_data = VtuReader.read_data

    def __init__(
        self,
        header_type,
        byte_order,
        compression,
        appended_data=None,
        raw_appended=None,
    ):
        if header_type not in ("UInt32", "UInt64"):
            raise ReadError(f"Unknown VTK XML header type '{header_type}'")
        if byte_order not in (None, "LittleEndian", "BigEndian"):
            raise ReadError(f"Unknown VTK XML byte order '{byte_order}'")
        if compression is not None:
            _compressor_for(compression)
        self.header_type = header_type
        self.byte_order = byte_order
        self.compression = compression
        self.appended_data = appended_data
        self.raw_appended = raw_appended


def load(filename, dataset_type, format_name):
    root, raw = _load_root(filename)
    if root.tag != "VTKFile":
        raise ReadError("Expected tag 'VTKFile'")
    if root.get("type") != dataset_type:
        raise ReadError(f"Expected type {dataset_type}")
    compression = root.get("compressor")
    # Preserve the structured readers' intentional lzma parity contract.
    # PolyData still has a Python-only lzma path, as it did before.
    if compression == "vtkLZMADataCompressor" and dataset_type != "PolyData":
        raise ReadError(f"lzma-compressed {format_name} is not supported")
    text = None
    appended = root.find("AppendedData")
    if appended is not None and raw is None:
        encoding = appended.get("encoding", "base64")
        if encoding != "base64":
            raise ReadError(f"Unknown {format_name} AppendedData encoding '{encoding}'")
        text = "".join((appended.text or "").split())
        if not text.startswith("_"):
            raise ReadError(f"{format_name}: AppendedData does not start with '_'")
        text = text[1:]
    reader = ArrayReader(
        root.get("header_type", "UInt32"),
        root.get("byte_order"),
        compression,
        text,
        raw,
    )
    return root, reader
