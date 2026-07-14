import { defineConfig } from "vitepress";

// https://vitepress.dev/reference/site-config
export default defineConfig({
  title: "meshio++",
  description: "I/O for many mesh formats",
  // Project site served from https://<org>.github.io/meshioplusplus/
  base: "/meshioplusplus/",
  lastUpdated: true,
  ignoreDeadLinks: true,

  themeConfig: {
    nav: [
      { text: "Quickstart", link: "/quickstart" },
      { text: "Formats", link: "/formats" },
      { text: "CLI", link: "/cli" },
    ],

    sidebar: [
      {
        text: "Introduction",
        items: [
          { text: "Overview", link: "/" },
          { text: "Installation", link: "/installation" },
          { text: "Quickstart", link: "/quickstart" },
        ],
      },
      {
        text: "Concepts",
        items: [
          { text: "Mesh data model", link: "/mesh_data_model" },
          { text: "Cell types", link: "/cell_types" },
        ],
      },
      {
        text: "Reference",
        items: [
          { text: "Supported formats", link: "/formats" },
          { text: "CLI reference", link: "/cli" },
          { text: "XDMF time series", link: "/xdmf_time_series" },
          { text: "Extending meshio++", link: "/extending" },
          { text: "ParaView plugin", link: "/paraview_plugin" },
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
          { text: "dolfin-xml", link: "/formats/dolfin" },
          { text: "exodus", link: "/formats/exodus" },
          { text: "flac3d", link: "/formats/flac3d" },
          { text: "flux", link: "/formats/flux" },
          { text: "freefem", link: "/formats/freefem" },
          { text: "gmsh", link: "/formats/gmsh" },
          { text: "h5m", link: "/formats/h5m" },
          { text: "hmf", link: "/formats/hmf" },
          { text: "mdpa", link: "/formats/mdpa" },
          { text: "med", link: "/formats/med" },
          { text: "medit", link: "/formats/medit" },
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
          { text: "ugrid", link: "/formats/ugrid" },
          { text: "unv", link: "/formats/unv" },
          { text: "vtk", link: "/formats/vtk" },
          { text: "vtu", link: "/formats/vtu" },
          { text: "wkt", link: "/formats/wkt" },
          { text: "xdmf", link: "/formats/xdmf" },
        ],
      },
    ],

    socialLinks: [
      { icon: "github", link: "https://github.com/<org>/meshioplusplus" },
    ],

    search: { provider: "local" },

    editLink: {
      pattern: "https://github.com/<org>/meshioplusplus/edit/main/doc/:path",
    },

    footer: {
      message: "Released under the MIT License.",
      copyright: "meshio++ contributors",
    },
  },
});
