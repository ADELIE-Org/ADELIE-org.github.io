# Packages

## Grid infrastructure

### [CartesianGrids.jl](https://github.com/ADELIE-org/CartesianGrids.jl)

[![docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://ADELIE-org.github.io/CartesianGrids.jl/dev)
[![CI](https://github.com/ADELIE-org/CartesianGrids.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/ADELIE-org/CartesianGrids.jl/actions)

Uniform Cartesian grid generation and management in arbitrary dimensions.

Provides node-based and cell-based grid abstractions (`CartesianGrid`, `SpaceTimeCartesianGrid`) backed by `StaticArrays` for type-stable coordinate arithmetic. A natural foundation for finite-difference and finite-volume solvers on structured meshes.

```julia
using CartesianGrids, StaticArrays

grid = CartesianGrid((-1.0, -1.0), (1.0, 1.0), (64, 64))
h    = meshsize(grid)          # node spacing
dx   = cell_spacing(grid)      # cell width
xc   = cell_center(grid, (3, 7))   # SVector coordinates
```

**Install:** `Pkg.add("CartesianGrids")`
