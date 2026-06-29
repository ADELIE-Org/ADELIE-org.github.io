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

#### [CartesianGeometry.jl](https://github.com/ADELIE-Org/CartesianGeometry.jl)

[![docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://ADELIE-Org.github.io/CartesianGeometry.jl/dev)
[![CI](https://github.com/ADELIE-Org/CartesianGeometry.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/ADELIE-Org/CartesianGeometry.jl/actions)

Cut-cell moment router — geometric moments of an implicit interface on a Cartesian grid, structured for the ADELIE operator algebra.

Given a `CartesianGrids` grid and an interface representation, produces the full moment set — volume fractions `V`, wet barycentres, face apertures `A`, second-kind staggered volumes `W` and apertures `B` — as plain `AbstractArray`s. Single dispatch on the interface type routes to either `:vofi` (libvofi C library, fast Float64) or `:vofijul` (pure-Julia, AD- and GPU-capable). N-dimensional, byte-identical to the Penguin router, up to ~30× faster and ~150× lower-allocation on the pure-Julia path.

```julia
using CartesianGeometry, CartesianGrids

disk(x, y) = sqrt(x^2 + y^2) - 0.5
grid = CartesianGrid((-1.0, -1.0), (1.0, 1.0), (65, 65))
ls   = LevelSetInterface(disk; backend = :vofi)

m = integrate(Tuple{0}, ls, grid, Float64, nan)   # volume moments
A = integrate(Tuple{1}, ls, grid, Float64, nan)   # face apertures
```

**Install:** `Pkg.add("CartesianGeometry")`

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

#### [ImplicitIntegration.jl](https://github.com/maltezfaria/ImplicitIntegration.jl) *(external)*

Numerical integration over implicitly-defined domains and surfaces — integrates functions over regions `φ(x) ≤ 0` and boundaries `φ(x) = 0` without requiring coordinate transformations.

**Install:** `Pkg.add(url = "https://github.com/maltezfaria/ImplicitIntegration.jl")`

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

Cut-cell geometric moments on a Cartesian grid from an explicit front — the front-tracking analogue of the level-set/VOF moment backends in the ADELIE geometry stack.

Takes a closed `CurveMesh` (2D) or `SurfaceMesh` (3D) from `FrontIntrinsicOps.jl` and returns, per cell, the wet volume fraction, barycentre, interface measure, cell type, face apertures, and second-kind staggered moments. Uses exact Sutherland–Hodgman polygon∩box clipping in 2D and its 3D analogue; converges O(h²) and matches the `vofi` backend element-by-element. Type-stable, allocation-free per cell, multithreaded, ForwardDiff-compatible for shape sensitivity.

```julia
using FrontCartesianGeometry, FrontIntrinsicOps

front = sample_circle(0.5, 512)
x = collect(range(-1.0, 1.0, length=129)); xyz = (x, copy(x))

V, bary, iface, celltype, _ = integrate(Tuple{0}, front, xyz, Float64, nan)
A = integrate(Tuple{1}, front, xyz, Float64, nan)   # face apertures
```

**Install:** `Pkg.add("FrontCartesianGeometry")`

#### [FrontIntrinsicOps.jl](https://github.com/ADELIE-Org/FrontIntrinsicOps.jl)

[![docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://ADELIE-Org.github.io/FrontIntrinsicOps.jl/dev)
[![CI](https://github.com/ADELIE-Org/FrontIntrinsicOps.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/ADELIE-Org/FrontIntrinsicOps.jl/actions)

Discrete exterior calculus (DEC) / FEEC toolkit on triangulated surface and polygonal curve meshes.

Provides geometry primitives (normals, areas, curvatures, Gauss–Bonnet), DEC/FEEC operators (incidence matrices `d₀`/`d₁`, Hodge stars, Laplace–Beltrami, Whitney forms, de Rham sequence), and surface PDE solvers (diffusion, advection–diffusion, heat-method geodesics, harmonic forms). Fully AD-compatible via ForwardDiff for shape sensitivity, multithreaded assembly, and a KernelAbstractions GPU backend.

```julia
using FrontIntrinsicOps

mesh = generate_icosphere(1.0, 3)
geom = compute_geometry(mesh)
dec  = build_dec(mesh, geom)

H = mean_curvature(mesh, geom, dec)   # ≈ 1/R on a unit sphere
u = solve_surface_helmholtz(mesh, geom, dec, f, 1.0)
```

**Install:** `Pkg.add("FrontIntrinsicOps")`

#### [FrontSpaceTimeDEC.jl](https://github.com/ADELIE-Org/FrontSpaceTimeDEC.jl)

[![docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://ADELIE-Org.github.io/FrontSpaceTimeDEC.jl/dev)
[![CI](https://github.com/ADELIE-Org/FrontSpaceTimeDEC.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/ADELIE-Org/FrontSpaceTimeDEC.jl/actions)

Space-time DEC for moving fronts — the time-dependent companion to `FrontIntrinsicOps.jl`.

Builds a one-slab worldsheet between two front states (`mesh_n` → `mesh_np1`): curve fronts sweep `(x,y,t) ∈ ℝ³`, surface fronts sweep `(x,y,z,t) ∈ ℝ⁴`. Provides space-time simplicial topology/incidence, barycentric dual measures and a diagonal Hodge star, and conservative moving-surface PDE assembly (mass, transport, diffusion, advection–diffusion) plus coupled nonlinear motion stepping. Type-stable, allocation-free on the hot geometric-measure path, ForwardDiff-differentiable through geometry → Hodge → assembly → solve (the linear solve routed through `FrontIntrinsicOps.linsolve`), multithreaded per-simplex fills, and a KernelAbstractions GPU path.

```julia
using FrontSpaceTimeDEC, FrontIntrinsicOps, StaticArrays

N = 64
ang   = [2π * (i - 1) / N for i in 1:N]
edges = [SVector(i, i == N ? 1 : i + 1) for i in 1:N]
mesh_n   = CurveMesh([SVector(cos(a), sin(a)) for a in ang], edges)
mesh_np1 = CurveMesh([SVector(1.02cos(a), 1.02sin(a)) for a in ang], edges)

slab = build_curve_slab(mesh_n, mesh_np1, 0.01)
geom = compute_spacetime_geometry(slab)
dec  = build_spacetime_dec(slab, geom)
u_np1, _ = solve_slab!(slab, geom, dec, [sin(2a) for a in ang]; diffusivity = 0.05)
```

**Install:** `Pkg.add("FrontSpaceTimeDEC")`
