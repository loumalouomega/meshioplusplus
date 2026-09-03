import { defineConfig } from "vitepress";

const SITE = "https://loumalouomega.github.io/meshioplusplus";
const DESCRIPTION =
  "Read and write 43 mesh formats through one data model, run 39 mesh and data operations on them, from Python, C, Fortran, Julia, R, WebAssembly or C++.";

// https://vitepress.dev/reference/site-config
export default defineConfig({
  title: "meshio++",
  description: "I/O and operations for many mesh formats",
  // Project site served from https://loumalouomega.github.io/meshioplusplus/
  base: "/meshioplusplus/",
  lastUpdated: true,
  ignoreDeadLinks: true,

  // Source-directory READMEs (TikZ icons, logo, diagram generator) are
  // developer notes, not site pages.
  srcExclude: ["icons/README.md", "logo/README.md", "diagrams/README.md"],

  head: [
    ["link", { rel: "icon", type: "image/svg+xml", href: "/meshioplusplus/logo-icon.svg" }],
    ["link", { rel: "alternate icon", type: "image/png", href: "/meshioplusplus/logo-icon.png" }],
    // Multi-resolution (16/32/48/64px) fallback for browsers/tools that only
    // look for a plain .ico (bookmark managers, old IE, some crawlers).
    ["link", { rel: "shortcut icon", type: "image/x-icon", href: "/meshioplusplus/favicon.ico" }],
    // Link previews (Slack, GitHub, Twitter/X, LinkedIn). The image is the
    // 1280x640 banner doc/logo/make_icon_assets.py generates; it is copied
    // into doc/public/ by hand like the other logo assets.
    ["meta", { property: "og:type", content: "website" }],
    ["meta", { property: "og:site_name", content: "meshio++" }],
    ["meta", { property: "og:title", content: "meshio++ — I/O and operations for many mesh formats" }],
    ["meta", { property: "og:description", content: DESCRIPTION }],
    ["meta", { property: "og:url", content: `${SITE}/` }],
    ["meta", { property: "og:image", content: `${SITE}/social-preview.png` }],
    ["meta", { property: "og:image:width", content: "1280" }],
    ["meta", { property: "og:image:height", content: "640" }],
    ["meta", { name: "twitter:card", content: "summary_large_image" }],
    ["meta", { name: "twitter:title", content: "meshio++" }],
    ["meta", { name: "twitter:description", content: DESCRIPTION }],
    ["meta", { name: "twitter:image", content: `${SITE}/social-preview.png` }],
  ],

  themeConfig: {
    logo: "/logo-icon.svg",

    nav: [
      { text: "Quickstart", link: "/quickstart" },
      { text: "Architecture", link: "/architecture" },
      { text: "Formats", link: "/formats" },
      { text: "CLI", link: "/cli" },
      { text: "Roadmap", link: "/roadmap" },
      { text: "Viewer", link: "/viewer/", target: "_self" },
      { text: "Dataset manager", link: "/viewer/dataset.html", target: "_self" },
      { text: "API (Doxygen)", link: "/api/", target: "_self" },
    ],

    sidebar: [
      {
        text: "Introduction",
        items: [
          { text: "Overview", link: "/" },
          { text: "Installation", link: "/installation" },
          { text: "Quickstart", link: "/quickstart" },
          { text: "Architecture", link: "/architecture" },
        ],
      },
      {
        text: "Concepts",
        items: [
          { text: "Mesh data model", link: "/mesh_data_model" },
          { text: "Cell types", link: "/cell_types" },
          { text: "Named regions", link: "/regions" },
          { text: "Polyhedra and ragged cells", link: "/polyhedra" },
          { text: "C++ mesh backends", link: "/cpp_backends" },
        ],
      },
      {
        text: "Reading and writing",
        items: [
          { text: "Supported formats", link: "/formats" },
          { text: "Selective reads", link: "/selective_read" },
          { text: "Memory-mapped reading", link: "/mmap" },
          { text: "Compression codecs", link: "/codecs" },
          { text: "XDMF time series", link: "/xdmf_time_series" },
          { text: "Provenance in written files", link: "/provenance" },
          { text: "Extending meshio++", link: "/extending" },
        ],
      },
      {
        text: "Mesh operations",
        items: [
          {
            text: "Inspection and topology",
            collapsed: true,
            items: [
              { text: "Mesh quality metrics", link: "/mesh_quality" },
              { text: "Geometric statistics", link: "/stats" },
              { text: "Mesh comparison (diff)", link: "/diff" },
              { text: "Surface extraction", link: "/extract_surface" },
              { text: "Skin extraction", link: "/extract_skin" },
              { text: "Reordering / renumbering", link: "/reorder" },
              { text: "Affine transform", link: "/transform" },
              { text: "Clean (weld / prune)", link: "/clean" },
              { text: "Merge / combine", link: "/merge" },
              { text: "Split (by criterion)", link: "/split" },
              { text: "Partitioning (N parts)", link: "/partition" },
              { text: "Cell conversion", link: "/convert_cells" },
            ],
          },
          {
            text: "Refinement and coarsening",
            collapsed: true,
            items: [
              { text: "Refinement (uniform and adaptive)", link: "/refine" },
              { text: "Green-element undo", link: "/undo_green" },
              { text: "Polyhedral refinement (subdivide)", link: "/subdivide" },
              { text: "Polyhedral coarsening (agglomerate)", link: "/agglomerate" },
              { text: "Decimation (QEM edge collapse)", link: "/decimate" },
              { text: "Volume decimation", link: "/decimate_volume" },
            ],
          },
          {
            text: "Remeshing and smoothing",
            collapsed: true,
            items: [
              { text: "Surface remeshing (ACVD)", link: "/remesh" },
              { text: "Volumetric remeshing (isosurface stuffing)", link: "/remesh_volume" },
              { text: "ODT remeshing (relocate + flip)", link: "/optimize_volume" },
              { text: "Smoothing (Laplacian / Taubin / ODT)", link: "/smooth" },
            ],
          },
          {
            text: "Fields and interpolation",
            collapsed: true,
            items: [
              { text: "Interpolation (field transfer)", link: "/interpolate" },
              { text: "Conservative interpolation", link: "/conservative_interpolate" },
              { text: "Field derivatives", link: "/gradient" },
              { text: "Second derivatives (Hessian)", link: "/hessian" },
              { text: "Error estimation", link: "/error" },
              { text: "Field integration", link: "/field_integration" },
            ],
          },
          {
            text: "Cutting, grids and distance",
            collapsed: true,
            items: [
              { text: "Crop (bbox / half-space / predicate)", link: "/crop" },
              { text: "Slicing / cross-sections", link: "/slice" },
              { text: "Isosurfaces / contours", link: "/isosurface" },
              { text: "Regular grids / voxelize", link: "/voxelize" },
              { text: "Signed distance", link: "/sdf" },
            ],
          },
        ],
      },
      {
        text: "Data operations",
        collapsed: true,
        items: [
          { text: "Overview", link: "/data_operations" },
          { text: "Array management", link: "/data_manage" },
          { text: "Location averaging", link: "/data_average" },
          { text: "Expressions (calc)", link: "/data_calc" },
          { text: "Value conditioning", link: "/data_condition" },
          { text: "Data summary", link: "/data_info" },
        ],
      },
      {
        text: "Pipelines and tools",
        items: [
          { text: "Settings pipeline", link: "/pipeline" },
          { text: "Sequences (transient)", link: "/sequences" },
          { text: "CLI reference", link: "/cli" },
          { text: "MCP server", link: "/mcp" },
          { text: "Interactive viewer", link: "/viewer" },
        ],
      },
      {
        text: "Integrations",
        items: [
          { text: "Interoperability", link: "/interop" },
          { text: "GPU handoff (DLPack / CuPy)", link: "/gpu" },
          { text: "ML data handling", link: "/ml" },
          { text: "Dataset manifests", link: "/datasets" },
          { text: "PhysicsNeMo integration", link: "/physicsnemo" },
          { text: "Blender add-on", link: "/blender" },
          { text: "ParaView plugin", link: "/paraview_plugin" },
        ],
      },
      {
        text: "Language surfaces",
        items: [
          { text: "WebAssembly / JavaScript", link: "/wasm" },
          { text: "C++ API", link: "/cpp_api" },
          { text: "C++ ABI compatibility", link: "/abi" },
          { text: "ABI additive-change reviews", link: "/abi_reviews" },
          { text: "Single-header C++", link: "/single_header" },
          { text: "C API", link: "/c_api" },
          { text: "Fortran", link: "/fortran" },
          { text: "Julia", link: "/julia" },
          { text: "R", link: "/r" },
          { text: "API reference (Doxygen)", link: "/api/", target: "_self" },
        ],
      },
      {
        text: "Project",
        items: [
          { text: "Benchmarks", link: "/benchmarks" },
          { text: "Roadmap", link: "/roadmap" },
        ],
      },
      {
        text: "Formats",
        collapsed: true,
        items: [
          { text: "abaqus", link: "/formats/abaqus" },
          { text: "ansys", link: "/formats/ansys" },
          { text: "ansysInp", link: "/formats/ansysinp" },
          { text: "avsucd", link: "/formats/avsucd" },
          { text: "cgns", link: "/formats/cgns" },
          { text: "dex", link: "/formats/dex" },
          { text: "dolfin-xml", link: "/formats/dolfin" },
          { text: "ensight", link: "/formats/ensight" },
          { text: "exodus", link: "/formats/exodus" },
          { text: "flac3d", link: "/formats/flac3d" },
          { text: "flux", link: "/formats/flux" },
          { text: "freefem", link: "/formats/freefem" },
          { text: "gid", link: "/formats/gid" },
          { text: "gmsh", link: "/formats/gmsh" },
          { text: "h5m", link: "/formats/h5m" },
          { text: "hmf", link: "/formats/hmf" },
          { text: "ip", link: "/formats/ip" },
          { text: "mdpa", link: "/formats/mdpa" },
          { text: "med", link: "/formats/med" },
          { text: "medit", link: "/formats/medit" },
          { text: "mff", link: "/formats/mff" },
          { text: "mfm", link: "/formats/mfm" },
          { text: "mphtxt", link: "/formats/mphtxt" },
          { text: "nastran", link: "/formats/nastran" },
          { text: "netgen", link: "/formats/netgen" },
          { text: "neuroglancer", link: "/formats/neuroglancer" },
          { text: "obj", link: "/formats/obj" },
          { text: "off", link: "/formats/off" },
          { text: "openfoam", link: "/formats/openfoam" },
          { text: "permas", link: "/formats/permas" },
          { text: "ply", link: "/formats/ply" },
          { text: "stl", link: "/formats/stl" },
          { text: "su2", link: "/formats/su2" },
          { text: "svg", link: "/formats/svg" },
          { text: "tecplot", link: "/formats/tecplot" },
          { text: "tetgen", link: "/formats/tetgen" },
          { text: "tikz", link: "/formats/tikz" },
          { text: "triangle", link: "/formats/triangle" },
          { text: "ugrid", link: "/formats/ugrid" },
          { text: "unv", link: "/formats/unv" },
          { text: "vti", link: "/formats/vti" },
          { text: "vtk", link: "/formats/vtk" },
          { text: "vtp", link: "/formats/vtp" },
          { text: "vtu", link: "/formats/vtu" },
          { text: "wkt", link: "/formats/wkt" },
          { text: "xdmf", link: "/formats/xdmf" },
        ],
      },
    ],

    socialLinks: [
      { icon: "github", link: "https://github.com/loumalouomega/meshioplusplus" },
    ],

    search: { provider: "local" },

    editLink: {
      // The repository's default branch is `master`; a `main` link 404s.
      pattern: "https://github.com/loumalouomega/meshioplusplus/edit/master/doc/:path",
    },

    footer: {
      message: "Released under the MIT License.",
      copyright: "meshio++ contributors",
    },
  },
});
