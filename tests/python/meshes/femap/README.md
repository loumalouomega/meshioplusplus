<!--pytest-codeblocks:skipfile-->
# Femap neutral reference fixtures

**Written by Femap, EMSolution and MYSTRAN, not by meshio++.** `tools/gen_femap_fixtures.py` assembles them from three MIT-licensed repositories (notices below). No Femap installation was available, so these files are the external-validation layer of the `femap` reader.

| File | Source | What it exercises |
|---|---|---|
| `A342.neu`, `A352.neu`, `A362.neu` | [FrontISTR](https://github.com/FrontISTR/FrontISTR) `fistr1/tools/neu2fstr/example/A/` (commit `dd7e2a08`), unmodified | Femap 8.2 tetra10, wedge15 and brick20: the degenerate-brick node slots |
| `A731.neu`, `B741.neu` | same, `A/` and `B/` | Femap 8.2 tri3 and quad4 |
| `MC361.neu` | same, `heat/` | Femap 8.2 brick8 with node groups (block 408) |
| `mesh_sample.neu` | [ems_file_format_converter](https://github.com/EMSolution-SSIL/ems_file_format_converter) `sample/`, unmodified | 4.41 file from EMSolution: tetra, pyramid (topology 14), wedge, brick |
| `ems_results_451.neu`, `ems_results_1051.neu` | the same repository: `mesh_sample.neu` followed by the result blocks of `post_sample451.neu` / `post_sample1051.neu` | 13 output sets written both as `451` records and as `1051` ranges |
| `v2020_results.neu`, `v82_results.neu`, `mystran_results.neu` | [femap_neutral_parser](https://framagit.org/numenic/femap_neutral_parser) `tests/data/` (commit `5951f191`): the header and result blocks of `FEMAP_v2020-1-0.neu`, `FEMAP_v8-2.neu` and `mystran_00.NEU`, with a 12-node, 11-bar model written by the generator between them | output sets and vectors written by Femap 2020.1 (`1051`, the seven-line output-set tail), Femap 8.2 (`451`) and MYSTRAN (`451`, a line before the first block) |

The results files carry no geometry of their own, so the generator adds the model they refer to: nodes 1-12 along x and bar `i` joining nodes `i` and `i+1`. Their values were checked against femap_neutral_parser's own reading of the originals (all 2,770 values of the 8.2 and MYSTRAN files equal); femap_neutral_parser cannot read the 2020.1 file.

## Licences

All three sources are MIT-licensed:

- FrontISTR: Copyright (c) 2019 FrontISTR Commons
- ems_file_format_converter: Copyright (c) 2025 Science Solutions International Laboratory, Inc.
- femap_neutral_parser: Copyright (c) 2021, Nicolas Cordier

```text
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```
