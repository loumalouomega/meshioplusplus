<!--pytest-codeblocks:skipfile-->
# LS-DYNA binout reference fixtures

| File | What it is |
|---|---|
| `binout_glstat` | a binout **LS-DYNA wrote** (`glstat` and `rwforc`, 101 outputs), copied unmodified from Ansys' [example data](https://github.com/ansys/example-data) (`result_files/binout_glstat`, MIT, notice below) |
| `nodout/binout` | the first 60 outputs of `nodout` and `elout/beam` and the first 7 of `glstat` and `matsum` of that repository's `result_files/binout_matsum` (a Hybrid III dummy, 15 MB), rewritten with lasso-python's LSDA writer by `tools/gen_binout_reference.py` |

`binout_reference.npz` freezes lasso-python 2.0.4's `Binout` reading of both (every numeric variable of every database), and `test_lsdyna_binout.py` compares both meshio++ engines against it; `tools/gen_binout_reference.py` regenerates the subset and the reference.

## Ansys example-data licence

MIT License

Copyright (c) 2020 - 2026 Synopsys, Inc. and ANSYS, Inc. All rights reserved.

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
