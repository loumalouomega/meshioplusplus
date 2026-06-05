# Installation

## Basic install

```
pip install meshio
```

The base install only requires NumPy. Most text-based formats work without any additional dependencies.

## Full install (all optional dependencies)

```
pip install meshio[all]
```

This pulls in:

| Package | Required for |
|---------|-------------|
| `h5py` | CGNS, H5M, HMF, MED, XDMF (HDF data format) |
| `netCDF4` | Exodus |

## Conda

```
conda install -c conda-forge meshio
```

## Development install

```
git clone https://github.com/nschloe/meshio.git
cd meshio
pip install -e ".[all]"
```

Run the test suite with:

```
pytest tests/
```

or via tox (tests against Python 3.8 and 3.12):

```
tox
```
