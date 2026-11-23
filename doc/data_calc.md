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

Functions: `abs(x)`, `sqrt(x)`, `min(a,b)`, `max(a,b)`, and `norm(v)` for the Euclidean magnitude of a vector, then the elementwise math functions:

| Function | Arguments | Result |
| --- | --- | --- |
| `exp`, `log`, `log10` | one | the same width as the argument |
| `sin`, `cos`, `tan`, `asin`, `acos`, `atan` | one (radians) | the same width |
| `sinh`, `cosh`, `tanh` | one | the same width |
| `floor`, `ceil`, `sign` | one | the same width (`sign` is −1, 0 or 1; NaN stays NaN) |
| `atan2(y, x)`, `pow(a, b)` | two | scalar-broadcasting, like `min` and `max` |

The **tuple functions** read a whole row and check its width, naming the function when it is wrong:

| Function | Operands | Result |
| --- | --- | --- |
| `dot(a, b)` | two of the same width | a scalar, `Σ aₖ·bₖ` |
| `cross(a, b)` | two 3-component | 3 components |
| `det(t)`, `trace(t)` | a 6-component symmetric `xx yy zz xy yz zx` or a 9-component row-major 3x3 | a scalar |
| `transpose(t)` | a 9-component 3x3 | 9 components |
| `mag(v)` | any | a scalar, the same as `norm(v)` |
| `matmul(a, b)` | a 9-component 3x3 and a 9-component 3x3 or a 3-component vector | 9 components (the product) or 3 (the matrix-vector product); the sum over `k` runs in ascending order |
| `outer(a, b)` | two 3-component | 9 components, `C[i][j] = a[i]·b[j]` |
| `inv(t)` | a 6-component symmetric or a 9-component 3x3 | the same width: `adj(t) / det(t)`, `NaN` in a row whose determinant is zero or non-finite |

Together with `pow` and `exp` they write most derived quantities in one call, for instance the helicity density `dot(U, curl)` from `gradient`'s `U:curl`, or `det(F)` of a deformation gradient, and `matmul(F, inv(F))` is the identity. `inv` is the same cofactor formula as [`tensor_invariants`](/tensor_invariants)' `inverse` output; there is no warning channel, so a singular row is `NaN` and nothing more. A 6-component tensor and a 9-component one are read as in [`tensor_invariants`](/tensor_invariants), and the determinant is the same triple product of rows.

Exponentiation is a **function**, `pow`, not an operator, so the grammar's precedence is unchanged. A domain error is a `NaN`, not an exception (`log(-1)`, `asin(2)`, `pow(-8, 0.5)`), and `log(0)` is `-inf`, the way division by zero is an infinity. There is no warning channel to count them in.

`abs`, `sqrt`, `floor`, `ceil`, `sign`, `min`, `max` and `norm` are exact, and the core and its Python twin agree on them bit for bit. `exp`, `log`, `log10`, the trigonometric and hyperbolic functions and `pow` come from the platform's libm in the core and from NumPy in the twin: they agree to about a relative `1e-12`, and may differ in the last bit between platforms, so do not compare their results for equality across engines or machines.

## Names

An identifier starts with a letter or `_` and may then contain letters, digits, `_`, `:` and `.`. The colon is load-bearing: it lets a mesh's own conventional names be referenced directly.

```python
mp.data_calc(mesh, "gmsh:physical + 1", location="cell", output="tag_plus_one")
```

A name containing spaces or operator characters can be written between backticks:

```python
mp.data_calc(mesh, "`my array` * 2", location="point", output="doubled")
```

### Coordinates

`x`, `y` and `z` (one component each) and `coords` (all three) name the geometry when no data array of that name exists:

```python
# A point field built from the coordinates, and from a cell's centroid:
mp.data_calc(mesh, "sqrt(x * x + y * y)", location="point", output="radius")
mp.data_calc(mesh, "z", location="cell", output="height")
```

- At `location="point"` they are the points' coordinates (a 2-D mesh has `z = 0`); at `location="cell"` they are the **centroid** of each cell, the mean of its nodes (the one [`partition`](/partition) and [`interpolate`](/interpolate) use), evaluated per cell block like every other cell operand.
- A **data array of that name wins**: an existing `x` array is read as an array, so no expression written before coordinates existed changes meaning. Backtick-quote nothing; rename the array or take the coordinate from `coords`.
- `location="field"` has no geometry to take them from: there `x` is an unknown array, as before.
- The numpy twin takes cell centroids from `partition`'s helper, which handles rectangular blocks only; a ragged or polyhedron block needs the C++ core.

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
unknown function 'foo' at position 0 (known: abs, sqrt, min, max, norm, exp, log, log10, sin, cos, tan, asin, acos, atan, atan2, sinh, cosh, tanh, floor, ceil, sign, pow, dot, cross, det, trace, transpose, mag, matmul, outer, inv)
'norm' takes exactly 1 argument (got 2) at position 0
'dot' needs two operands of the same width (got 3 and 9) at position 0
'det' needs a 6-component (symmetric) or 9-component (3x3) operand (got 3) at position 0
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
