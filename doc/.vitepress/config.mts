import { defineConfig } from "vitepress";

// https://vitepress.dev/reference/site-config
export default defineConfig({
  title: "meshio",
  description: "I/O for many mesh formats",
  // Project site served from https://<user>.github.io/meshio/
  base: "/meshio/",
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
          { text: "Extending meshio", link: "/extending" },
          { text: "ParaView plugin", link: "/paraview_plugin" },
        ],
      },
    ],

    socialLinks: [
      { icon: "github", link: "https://github.com/nschloe/meshio" },
    ],

    search: { provider: "local" },

    editLink: {
      pattern: "https://github.com/nschloe/meshio/edit/main/doc/:path",
    },

    footer: {
      message: "Released under the MIT License.",
      copyright: "meshio contributors",
    },
  },
});
