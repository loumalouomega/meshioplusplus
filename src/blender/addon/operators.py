"""The File > Import and File > Export operators.

Thin by design: every decision about what a mesh becomes lives in
``meshioplusplus._blender``, which ships in the wheel and is tested without
Blender. These classes gather options, call it, and translate a failure into a
report the user can read.
"""

# NOTE: no `from __future__ import annotations` here, and that is load-bearing.
# Blender registers properties by reading `__annotations__` at class creation
# and needs the real `_PropertyDeferred` object; under PEP 563 every one of
# them becomes a string and registration silently produces an operator with no
# properties at all.

import os

import bpy
from bpy.props import (
    BoolProperty,
    CollectionProperty,
    EnumProperty,
    FloatProperty,
    StringProperty,
)
from bpy_extras.io_utils import (
    ExportHelper,
    ImportHelper,
    axis_conversion,
    orientation_helper,
)
from mathutils import Matrix

from . import formats


def _report_errors(method):
    """Turn any exception into an operator error report.

    ``_server.py``'s rule -- a tool never surfaces a traceback -- applied to a
    UI. A malformed file is an ordinary event here, not a crash.
    """

    def wrapper(self, context):
        try:
            return method(self, context)
        except Exception as exc:  # noqa: BLE001 - the UI boundary
            self.report({"ERROR"}, f"{type(exc).__name__}: {exc}")
            return {"CANCELLED"}

    wrapper.__name__ = method.__name__
    wrapper.__doc__ = method.__doc__
    return wrapper


@orientation_helper(axis_forward="Y", axis_up="Z")
class MESHIOPLUSPLUS_OT_import_mesh(bpy.types.Operator, ImportHelper):
    """Import a mesh through meshio++"""

    bl_idname = "meshioplusplus.import_mesh"
    bl_label = "Import mesh (meshio++)"
    bl_options = {"REGISTER", "UNDO", "PRESET"}

    filename_ext = ""
    filter_glob: StringProperty(default=formats.import_glob(), options={"HIDDEN"})
    files: CollectionProperty(
        type=bpy.types.OperatorFileListElement, options={"HIDDEN", "SKIP_SAVE"}
    )
    directory: StringProperty(subtype="DIR_PATH", options={"HIDDEN", "SKIP_SAVE"})

    show_all_files: BoolProperty(
        name="Show all files",
        description=(
            "Ignore the extension filter. meshio++ knows more extensions than "
            "Blender's filter field can hold, so a few are not listed by default"
        ),
        default=False,
    )
    file_format: EnumProperty(
        name="Format",
        items=lambda self, context: formats.read_items(),
        description="Force a reader instead of inferring one from the file name",
    )
    global_scale: FloatProperty(
        name="Scale", default=1.0, min=1e-6, max=1e6, soft_min=0.001, soft_max=1000.0
    )
    geometry_only: BoolProperty(
        name="Geometry only",
        description="Skip every data array (much faster for a quick look)",
        default=False,
    )
    import_point_data: BoolProperty(name="Point data", default=True)
    import_cell_data: BoolProperty(name="Cell data", default=True)
    import_regions: BoolProperty(name="Regions", default=True)
    import_field_data: BoolProperty(name="Field data", default=True)
    validate_mesh: BoolProperty(
        name="Validate",
        description="Let Blender check and repair the imported topology",
        default=True,
    )

    def draw(self, context):
        layout = self.layout
        layout.use_property_split = True
        layout.prop(self, "file_format")
        layout.prop(self, "global_scale")
        layout.prop(self, "show_all_files")
        column = layout.column(heading="Import")
        column.prop(self, "geometry_only")
        sub = column.column()
        sub.enabled = not self.geometry_only
        sub.prop(self, "import_point_data")
        sub.prop(self, "import_cell_data")
        sub.prop(self, "import_regions")
        sub.prop(self, "import_field_data")
        layout.prop(self, "validate_mesh")

    def invoke(self, context, event):
        if self.show_all_files:
            self.filter_glob = ""
        return ImportHelper.invoke(self, context, event)

    def _paths(self):
        if self.files and self.directory:
            return [os.path.join(self.directory, f.name) for f in self.files if f.name]
        return [self.filepath]

    def _matrix(self):
        # Axis conversion and scale ride the object matrix, never the vertex
        # array: free, non-destructive, undoable, and it keeps the conversion
        # layer free of transform code (meshio++ has `transform` for baking).
        rotation = axis_conversion(
            from_forward=self.axis_forward, from_up=self.axis_up
        ).to_4x4()
        return rotation @ Matrix.Scale(self.global_scale, 4)

    @_report_errors
    def execute(self, context):
        import meshioplusplus as mio

        matrix = self._matrix()
        imported = 0
        for path in self._paths():
            mesh = mio.read(
                path,
                file_format=None if self.file_format == "AUTO" else self.file_format,
                arrays=[] if self.geometry_only else None,
            )
            data = mio.to_blender(
                mesh,
                name=os.path.splitext(os.path.basename(path))[0],
                point_data=self.import_point_data and not self.geometry_only,
                cell_data=self.import_cell_data and not self.geometry_only,
                regions=self.import_regions and not self.geometry_only,
                field_data=self.import_field_data and not self.geometry_only,
                validate=self.validate_mesh,
            )
            obj = bpy.data.objects.new(data.name, data)
            obj.matrix_world = matrix
            context.collection.objects.link(obj)
            imported += 1
        self.report({"INFO"}, f"meshio++: imported {imported} mesh(es)")
        return {"FINISHED"}


@orientation_helper(axis_forward="Y", axis_up="Z")
class MESHIOPLUSPLUS_OT_export_mesh(bpy.types.Operator, ExportHelper):
    """Export the scene through meshio++"""

    bl_idname = "meshioplusplus.export_mesh"
    bl_label = "Export mesh (meshio++)"
    bl_options = {"REGISTER", "PRESET"}

    filename_ext = ".vtu"
    filter_glob: StringProperty(default=formats.import_glob(), options={"HIDDEN"})

    file_format: EnumProperty(
        name="Format",
        items=lambda self, context: formats.write_items(),
        description="Force a writer instead of inferring one from the file name",
    )
    use_selection: BoolProperty(name="Selection only", default=False)
    apply_modifiers: BoolProperty(name="Apply modifiers", default=True)
    export_attributes: BoolProperty(name="Attributes", default=True)
    global_scale: FloatProperty(name="Scale", default=1.0, min=1e-6, max=1e6)

    def draw(self, context):
        layout = self.layout
        layout.use_property_split = True
        layout.prop(self, "file_format")
        layout.prop(self, "use_selection")
        layout.prop(self, "apply_modifiers")
        layout.prop(self, "export_attributes")
        layout.prop(self, "global_scale")

    @_report_errors
    def execute(self, context):
        import meshioplusplus as mio

        source = (
            context.selected_objects if self.use_selection else context.scene.objects
        )
        objects = [obj for obj in source if obj.type == "MESH"]
        if not objects:
            self.report({"ERROR"}, "meshio++: no mesh object to export")
            return {"CANCELLED"}

        conversion = axis_conversion(
            to_forward=self.axis_forward, to_up=self.axis_up
        ).to_4x4() @ Matrix.Scale(self.global_scale, 4)

        meshes = []
        for obj in objects:
            mesh = mio.from_blender(
                obj,
                apply_modifiers=self.apply_modifiers,
                attributes=self.export_attributes,
            )
            # Bake each object's own placement, which a file format has no way
            # to carry separately from the coordinates.
            matrix = conversion @ obj.matrix_world
            meshes.append(mio.transform(mesh, matrix=[list(row) for row in matrix]))

        # More than one object is `merge`, not a hand-rolled concatenation.
        mesh = meshes[0] if len(meshes) == 1 else mio.merge(meshes)
        mio.write(
            self.filepath,
            mesh,
            file_format=None if self.file_format == "AUTO" else self.file_format,
        )
        self.report({"INFO"}, f"meshio++: exported {len(objects)} object(s)")
        return {"FINISHED"}


CLASSES = (MESHIOPLUSPLUS_OT_import_mesh, MESHIOPLUSPLUS_OT_export_mesh)


def menu_import(self, context):
    self.layout.operator(
        MESHIOPLUSPLUS_OT_import_mesh.bl_idname, text="Mesh via meshio++"
    )


def menu_export(self, context):
    self.layout.operator(
        MESHIOPLUSPLUS_OT_export_mesh.bl_idname, text="Mesh via meshio++"
    )
