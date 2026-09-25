# Tecplot `.szplt`

Three ways, from none to all of Tecplot:

1. **Read it directly** with a meshio++ built against TecIO, Tecplot's freely available I/O library: the [`szplt`](../formats/szplt.md) format. No Tecplot licence is needed.
2. **PyTecplot**, with a licensed Tecplot 360: `contrib/tecplot/szplt_to_plt.py` loads each `.szplt` (`tecplot.data.load_tecplot_szl`) and saves it as binary `.plt` (`tecplot.data.save_tecplot_plt`), which meshio++ reads natively ([Tecplot](../formats/tecplot.md)):

   ```bash
   python szplt_to_plt.py run.szplt more.szplt           # run.plt, more.plt
   ```

3. **A Tecplot macro**, `contrib/tecplot/szplt_to_plt.mcr`: edit its two paths, then *Scripting → Play Macro/Script*, or `tec360 -b -p szplt_to_plt.mcr`. It loads the file with the Subzone Data Loader and writes it with `$!WRITEDATASET … BINARY = YES`, keeping the data-sharing links.

The first is tested on TecIO-written files; the PyTecplot script against a stand-in (`tests/python/test_contrib_routes.py`); the run against Tecplot 360 itself — the macro included — is listed in the [roadmap](../roadmap.md#awaiting-a-licensed-run).
