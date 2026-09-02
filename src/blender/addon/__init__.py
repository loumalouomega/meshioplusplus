"""meshio++ for Blender: import and export 40+ FEA and scientific mesh formats.

Deliberately no ``bl_info``: this is a Blender 4.2+ *extension*, so
``blender_manifest.toml`` carries the metadata and a leftover ``bl_info`` is a
validation warning.

All of the conversion lives in ``meshioplusplus._blender``, which arrives as a
bundled wheel. This package is only the Blender-facing shell.
"""

import bpy

from . import operators


def register():
    for cls in operators.CLASSES:
        bpy.utils.register_class(cls)
    bpy.types.TOPBAR_MT_file_import.append(operators.menu_import)
    bpy.types.TOPBAR_MT_file_export.append(operators.menu_export)


def unregister():
    bpy.types.TOPBAR_MT_file_export.remove(operators.menu_export)
    bpy.types.TOPBAR_MT_file_import.remove(operators.menu_import)
    for cls in reversed(operators.CLASSES):
        bpy.utils.unregister_class(cls)
