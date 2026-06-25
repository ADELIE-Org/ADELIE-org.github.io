# Packages

## Grid infrastructure

#### [CartesianGrids.jl](https://github.com/ADELIE-org/CartesianGrids.jl)

[![docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://ADELIE-org.github.io/CartesianGrids.jl/dev)
[![CI](https://github.com/ADELIE-org/CartesianGrids.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/ADELIE-org/CartesianGrids.jl/actions)

Uniform Cartesian grid generation and management in arbitrary dimensions.

Provides node-based and cell-based grid abstractions (`CartesianGrid`, `SpaceTimeCartesianGrid`) backed by `StaticArrays` for type-stable coordinate arithmetic. A natural foundation for finite-difference and finite-volume solvers on structured meshes.

```julia
using CartesianGrids, StaticArrays

grid = CartesianGrid((-1.0, -1.0), (1.0, 1.0), (64, 64))
h    = meshsize(grid)              # node spacing
dx   = cell_spacing(grid)          # cell width
xc   = cell_center(grid, (3, 7))   # SVector coordinates
```

**Install:** `Pkg.add("CartesianGrids")`

## Geometry

### Signed Distance / Volume-of-Fluid Initialization

#### [Vofinit.jl](https://github.com/ADELIE-org/Vofinit.jl)

[![docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://ADELIE-org.github.io/Vofinit.jl/dev)
[![CI](https://github.com/ADELIE-org/Vofinit.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/ADELIE-org/Vofinit.jl/actions)

Volume-of-fluid (VOF) initialization — exact cell types and volume fractions from an implicit interface.

Julia wrapper around the [VOFI](https://github.com/VOFTracking/vofi) library (`libvofi_jll`) for computing exact cell volume and type from a level-set function. Provides allocation-free callbacks for use in hot loops.

```julia
using Vofinit

x0 = Cdouble.((0.0, 0.0, 0.0))
h0 = Cdouble.((1.0, 1.0, 1.0))

celltype = getcelltype(x0, h0) do x...
    sum(x .* x) - 0.25          # implicit sphere of radius 0.5
end

vf = getcc(x0, h0) do x...
    sum(x .* x) - 0.25
end
```

**Install:** `Pkg.add("Vofinit")`

#### [VofiJul.jl](https://github.com/ADELIE-Org/VofiJul.jl)

[![docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://ADELIE-Org.github.io/VofiJul.jl/dev)
[![CI](https://github.com/ADELIE-Org/VofiJul.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/ADELIE-Org/VofiJul.jl/actions)

Pure-Julia VOFI — volume fractions, interface moments, and adjoints in 1–4D.

Full port of the VOFI algorithm: `vofi_get_cc` (fraction + centroid + interface measure), `vofi_get_cell_type`, and `vofi_interface_centroid`. AD-compatible via ForwardDiff for geometry sensitivity and a custom shape-derivative adjoint for level-set parameter sensitivity. Zero-allocation hot paths via a per-thread cached workspace; thread-safe batch API and a KernelAbstractions GPU kernel.

```julia
using VofiJul

vf = vofi_get_cc((x, y, z) -> x^2 + y^2 + z^2 - 0.25,
                 (0.0, 0.0, 0.0), (1.0, 1.0, 1.0))
```

#### [ImplicitIntegration.jl](https://github.com/ELallinec/ImplicitIntegration.jl) *(external)*

Numerical integration over implicitly-defined domains and surfaces — integrates functions over regions `φ(x) ≤ 0` and boundaries `φ(x) = 0` without requiring coordinate transformations.

**Install:** `Pkg.add(url = "https://github.com/ELallinec/ImplicitIntegration.jl")`

### VOF

**Install:** `Pkg.dev("VofiJul")`

#### [VOFTools.jl](https://github.com/ADELIE-Org/VOFTools.jl)

[![docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://ADELIE-Org.github.io/VOFTools.jl/dev)
[![CI](https://github.com/ADELIE-Org/VOFTools.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/ADELIE-Org/VOFTools.jl/actions)

Geometrical tools for VOF methods: polygon/polyhedron truncation, volume-fraction computation, and PLIC reconstruction.

Pure-Julia port of [VOFTools 5](https://www.dimf.upct.es/personal/lopez_jp/VOFTools/) (López & Hernández). Provides convex-cell truncation (`inte2d!`, `inte3d!`), area/volume computation, volume enforcement for general and rectangular cells (CIBRAVE), and material initialisation from an implicit interface. Zero-allocation hot paths, thread-safe, and AD-compatible.

```julia
using VOFTools

# Truncate a unit square by the half-plane x ≤ 0.6 and get the cut area
verts = [0.0 1.0 1.0 0.0; 0.0 0.0 1.0 1.0]
n, area = inte2d!(verts, [1.0, 0.0], 0.6)
```

**Install:** `Pkg.add(url = "https://github.com/ADELIE-Org/VOFTools.jl")`

#### [isoap.jl](https://github.com/ADELIE-Org/isoap.jl)

[![docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://ADELIE-Org.github.io/isoap.jl/dev)
[![CI](https://github.com/ADELIE-Org/isoap.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/ADELIE-Org/isoap.jl/actions)

Isosurface extraction on arbitrary polyhedra and grids — pure-Julia port of the ISOAP algorithm (López & Hernández, 2020).

Reconstructs the surface `φ(x) = φ_iso` inside general polyhedral cells. Supports 15 built-in cell geometries, structured and OpenFOAM grids with VTK output, differentiable coordinates (ForwardDiff-compatible), a multithreaded grid extractor, and a KernelAbstractions kernel for GPU execution.

```julia
using ISOAP

ws  = IsoapWorkspace()
pol = isoap!(ws, cube(), φ_values, 0.0)   # extract iso-polygon at φ = 0
```

**Install:** `Pkg.dev("ISOAP")`

### Lagrangian tracking

#### [FrontCartesianGeometry.jl](https://github.com/ADELIE-Org/FrontCartesianGeometry.jl)

[![docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://ADELIE-Org.github.io/FrontCartesianGeometry.jl/dev)
[![CI](https://github.com/ADELIE-Org/FrontCartesianGeometry.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/ADELIE-Org/FrontCartesianGeometry.jl/actions)

Cartesian front-tracking geometry primitives.

**Install:** `Pkg.dev("FrontCartesianGeometry")`

#### [FrontIntrinsicOps.jl](https://github.com/ADELIE-Org/FrontIntrinsicOps.jl)

[![docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://ADELIE-Org.github.io/FrontIntrinsicOps.jl/dev)
[![CI](https://github.com/ADELIE-Org/FrontIntrinsicOps.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/ADELIE-Org/FrontIntrinsicOps.jl/actions)

Intrinsic differential operators on fronts and interfaces.

**Install:** `Pkg.dev("FrontIntrinsicOps")`
