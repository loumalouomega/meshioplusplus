# Data expressions (`data_calc`)

`data_calc` derives a new data array from a small elementwise expression over the existing arrays at the same location. It is a [data operation](/data_operations), not a file format; the geometry is never modified.

![A derived speed field, `norm(velocity)`, coloured on the mesh](/images/data_calc_speed.png)

```python
import meshioplusplus as mp

mesh = mp.read("flow.vtu")

mesh = mp.data_calc(mesh, "norm(velocity)", location="point", output="speed")
mesh = mp.data_calc(mesh, "p_new - p_old", location="cell", output="dp")
mesh = mp.data_calc(mesh, "0.5 * rho * speed * speed", location="point",
                    output="dynamic_pressure")
```

::: info No arbitrary code
The evaluator is a hand-written tokenizer plus recursive-descent parser restricted to the grammar below. There is **no external parser library** and **no evaluation of arbitrary code** — an expression that is not built from these tokens is a diagnosed error, never something that executes.
:::

## Grammar

```
expr    := term (('+'|'-') term)*
term    := unary (('*'|'/') unary)*
unary   := ('-'|'+') unary | primary
primary := number | ident | ident '(' expr (',' expr)* ')' | '(' expr ')'
```

Functions: `abs(x)`, `sqrt(x)`, `min(a,b)`, `max(a,b)`, and `norm(v)` for the Euclidean magnitude of a vector. Exponentiation is deliberately absent.

## Names

An identifier starts with a letter or `_` and may then contain letters, digits, `_`, `:` and `.`. The colon is load-bearing: it lets a mesh's own conventional names be referenced directly.

```python
mp.data_calc(mesh, "gmsh:physical + 1", location="cell", output="tag_plus_one")
```

A name containing spaces or operator characters can be written between backticks:

```python
mp.data_calc(mesh, "`my array` * 2", location="point", output="doubled")
```

## Values and broadcasting

Every operand is either a scalar (1 component) or a vector/tensor of up to 16 components. Binary operators require matching component counts, except that a **scalar broadcasts** against a wider operand:

```python
mp.data_calc(mesh, "velocity * 2",  ...)   # (n, 3) -> (n, 3)
mp.data_calc(mesh, "norm(velocity)", ...)  # (n, 3) -> (n,)   norm is always scalar
mp.data_calc(mesh, "velocity - velocity", ...)  # (n, 3) -> (n, 3)
```

Combining, say, a 3-component with a 9-component array is an error. A result with one component is stored 1-D; otherwise it keeps its `(n, ncomp)` shape.

Arithmetic is always performed in `double` and the result is stored as `float64`. Division by zero yields an IEEE infinity or `NaN` and is **not** an error — see the [NaN policy](/data_operations#nan-policy).

For `location="cell"` the expression is evaluated **once per cell block**, so the produced `cell_data` always has exactly one array per block.

## Errors

Every failure raises `ValueError`, prefixed `meshio++: data_calc: ` and, where meaningful, carrying the 0-based character position within the expression:

```
unexpected character '#' at position 7
unexpected end of expression
expected ')' at position 12
trailing input after the expression at position 9
unknown point_data array 'temp' at position 5 (available: T, p, u)
unknown function 'log' at position 0 (known: abs, sqrt, min, max, norm)
'norm' takes exactly 1 argument (got 2) at position 0
'min' takes exactly 2 arguments (got 1) at position 3
cannot combine a 3-component array with a 9-component array in '*' at position 12
array 'u' has 42 rows but point_data needs 100
expression nests deeper than 64 levels
result has 24 components, more than the supported maximum of 16
output name 'T' already exists in point_data (pass overwrite=true to replace it)
```

The nesting limit guards the recursive-descent parser against a stack overflow from hostile input, since expressions reach this code from the C ABI and both CLIs.

## CLI

```bash
meshioplusplus data calc in.vtu out.vtu --point "speed = norm(velocity)"
meshioplusplus data calc in.vtu out.vtu --cell  "dp = p_new - p_old"
```

The value is split on the **first** `=` into the output name and the expression. `--point`, `--cell` and `--field` are repeatable, so several arrays can be derived in one invocation; pass `--overwrite` to allow replacing an existing array. See the [CLI reference](/cli#meshioplusplus-data).

## Other languages

- **C API** — `mio_data_calc(mesh, expression, location, output_name, overwrite)` returns a new `mio_mesh*`, or `NULL` with `mio_last_error()` set. See the [C API reference](/c_api).
- **Fortran** — `m%data_calc('norm(velocity)', 'speed')`, with optional `location=` and `overwrite=`. See the [Fortran reference](/fortran).
- **WebAssembly / JavaScript** — `dataCalc(mesh, "norm(velocity)", "point", "speed", false)`; a bad expression throws a catchable `Error`. See the [WebAssembly reference](/wasm).
