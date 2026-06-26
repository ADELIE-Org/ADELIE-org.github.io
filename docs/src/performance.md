# Performance & gains vs. the original Penguin packages

The ADELIE geometry packages are oracle-guarded ports of the original
[Penguin](https://github.com/Fastaxx) cut-cell codes (and their Fortran/C
ancestors). Beyond reproducing the original results bit-for-bit, the ports were
refactored along five axes — **type-stability, allocation, multithreading,
automatic differentiation (AD), and GPU**. This page records the measured gains.

!!! note "How to read these numbers"
    *Allocation* figures are **measured per cut cell** (warmed, scratch reused),
    using Julia's `@allocated`. The "original" column is the Penguin baseline
    captured by importing the unmodified code before refactoring. *Runtime* gains
    follow directly from removing per-cell heap traffic (garbage-collector
    pressure) and were measured on the same cells / test suites. *Threading* is
    reported as validated determinism plus observed scaling — each package is
    reentrant and bit-for-bit identical serial vs. multi-threaded.

## Allocation — the headline

[`VofiJul`](https://github.com/ADELIE-Org/VofiJul.jl) (the Julia VOFI level-set
cut-cell engine) carried by far the largest reduction: its original per-cut-cell
allocation grew steeply with dimension. After the refactor (reused per-thread
workspace, concrete integrand wrapper, stack-resident `SArray` stencils,
allocation-free static-array resets), per-cut-cell allocation is **near the
floor in every dimension**.

| Path (per cut cell) | Original (Penguin) | ADELIE | Reduction |
|---|--:|--:|--:|
| **2D** volume fraction | ≈ 81 KB | **208 B** | ≈ 390× |
| **3D** volume fraction | ≈ 2.59 MB | **304 B** | ≈ 8500× |
| **4D** hypervolume fraction | ≈ 48 MB | **≈ 6 KB** | ≈ 8000× |
| Interface centroid (2D/3D) | ≈ 6.2 KB | **0 B** | ∞ |

Full / empty (non-cut) cells allocate ≈ 0 in all dimensions, and a plain
`for cell; vofi_get_cc(...)` loop now allocates ≈ 0 per cell (not just the
threaded batch driver), because the workspace is fetched from a per-thread cache.

[`Vofinit`](https://github.com/ADELIE-Org/Vofinit.jl) (the `libvofi` C-library
binding) had a different but equally important reduction. Its integrand callback
used to take an untyped `Any` thunk, so the `Cdouble` result was **boxed once per
quadrature-point evaluation** — allocation *grew with the quadrature order* `npt`
(≈ 13 KB+ per cut cell). Routing the callback through a concretely-typed
`FunctionWrapper` (`VofiCB`) makes the per-evaluation call allocation-free:

| Path (per cut cell) | Original | ADELIE | Property |
|---|--:|--:|---|
| `Vofinit` `getcc` | ≈ 13 KB+ (∝ `npt`) | **≈ 80 B** | now **constant**, independent of quadrature order |

### The other geometry backends were already lean

The remaining two backends did **not** need an allocation rescue — their ports
were allocation-free by design — so their gain is the *new capability* (AD /
threading / GPU) rather than a reduction:

| Package | 2D | 3D | 4D | Note |
|---|--:|--:|--:|---|
| [`isoap`](https://github.com/ADELIE-Org/isoap.jl) | 0 B | 0 B | — | `isoap!`/`isopol!` 0-alloc; the GPU kernel core `isoap_cell!` is also 0-alloc (stack-only, device-ready) |
| [`VOFTools`](https://github.com/ADELIE-Org/VOFTools.jl) | 0 B | 0 B | — | `inte2d!`/`inte3d!`/`enforv*` 0-alloc (TLS-pooled scratch) |

## CartesianGeometry router — end-to-end moments

[`CartesianGeometry`](https://github.com/ADELIE-Org/CartesianGeometry.jl) is the
moment-producing **router**: it loops over the grid, classifies each cell
(full / empty / cut), dispatches cut cells to a backend, and assembles the
moment arrays. The numbers below are *whole-router* `integrate(...)` calls (loop +
classification + executor + container assembly), not isolated per-cell kernels —
the end-to-end cost a caller actually pays. Measured old Penguin
`CartesianGeometry` vs. new ADELIE on the **same level set** (sphere, radius 0.5)
and grids — **2D = 129², 3D = 33³** — best-of-7 wall time, warmed `@allocated`.
The four moment families: `V`/centroid/interface come from one `Tuple{0}` call;
`A` = first-kind face apertures (`Tuple{1}`); `W` = second-kind staggered volumes
(`Tuple{0}` + barycenters); `B` = second-kind staggered apertures (`Tuple{1}` +
barycenters).

### `:vofijul` backend (pure-Julia executor — the refactor payoff)

This is the path that was rebuilt for type-stability / allocation / threading /
AD / GPU. Old-vs-new:

| Dim | Moment | Old time | New time | Speedup | Old alloc | New alloc | Reduction |
|---|---|--:|--:|--:|--:|--:|--:|
| 2D | `V`/centroid/interface | 54.4 ms | **1.0 ms** | ≈ 52× | 50.4 MB | **0.94 MB** | ≈ 54× |
| 2D | `A` | 0.93 ms | **0.15 ms** | ≈ 6× | 1.87 MB | **0.27 MB** | ≈ 7.3× |
| 2D | `W` | 26.8 ms | **1.99 ms** | ≈ 13× | 25.5 MB | **0.39 MB** | ≈ 69× |
| 2D | `B` | 1.11 ms | **0.16 ms** | ≈ 7× | 2.36 MB | **0.27 MB** | ≈ 9.3× |
| 3D | `V`/centroid/interface | 545 ms | **14.8 ms** | ≈ 37× | 405 MB | **2.79 MB** | ≈ 145× |
| 3D | `A` | 109 ms | **4.9 ms** | ≈ 22× | 65.0 MB | **1.34 MB** | ≈ 51× |
| 3D | `W` | 1423 ms | **50.6 ms** | ≈ 28× | 1038 MB | **2.03 MB** | ≈ **537×** |
| 3D | `B` | 133 ms | **5.6 ms** | ≈ 24× | 78.9 MB | **1.46 MB** | ≈ 57× |

The volume-type moments `V`/`W` carry the full VofiJul executor (up to
**537× less allocation in 3D**); the apertures `A`/`B` were freed of a
per-cut-cell boxing cost in the router engine (see below), so they too dropped
by an order of magnitude in allocation and ran several times faster.

### `:vofi` backend (libvofi C executor)

The C executor itself is unchanged, so every gain here is *router-level*:
cut-only interface-centroid evaluation, a reflection-free integrand wrapper,
reused buffers, and — the largest lever — **specializing the moment engine on
the level-set function** (`f::F` type parameter on the `Tuple{1}` / second-kind
`integrate!` methods, matching what the `@generated` `Tuple{0}` path already
did). Without it the engine declined to specialize on the `Function` argument,
heap-boxing the per-cell `SVector` stencils:

| Dim | Moment | Speedup | Alloc reduction |
|---|---|--:|--:|
| 2D | `V`/centroid/interface | ≈ 1.1× | ≈ 1.8× |
| 2D | `A` | ≈ 5× | ≈ 7.3× |
| 2D | `W` | ≈ 1.4× | ≈ 16× |
| 2D | `B` | ≈ 6× | ≈ 9.3× |
| 3D | `V`/centroid/interface | ≈ 1.5× | ≈ 9.3× |
| 3D | `A` | ≈ 2.1× | ≈ 14× |
| 3D | `W` | ≈ 1.8× | ≈ 79× |
| 3D | `B` | ≈ 2.0× | ≈ 15× |

Every moment remains **byte-identical** to old Penguin (`max|Δ| = 0`
element-by-element across all fields, both backends); these are pure performance
gains, not a change in output. The `:implicitintegration` backend got the same
treatment on its face path (promoting the fixed face coordinate to the
quadrature element type to avoid a `Union`-typed integrand argument), cutting its
aperture allocation by ≈ 70%.

## FrontCartesianGeometry — Lagrangian front → cut-cell moments

[`FrontCartesianGeometry`](https://github.com/ADELIE-Org/FrontCartesianGeometry.jl)
is the **front-driven** executor: instead of a signed-distance level set it takes
an explicit Lagrangian marker mesh (a polygon in 2D, a triangulated surface in
3D) and computes the same moment families by *clipping* each cut cell against the
front. It is the path you use when you only have markers, not an analytic φ.

The five-axis refactor (per-cut-cell `FaceLoopWorkspace`, 2D lazy-clip/ping-pong,
3D Sutherland–Hodgman box-face clip with CSR-flattened polyhedron loops) produced
the largest allocation collapse of any package. Old Penguin
`FrontCartesianGeometry` vs. new ADELIE on the **same fronts** —
`sample_circle(0.5, 512)` (2D) and `generate_icosphere(0.5, 2)` (3D) — same grids
(**2D = 129², 3D = 33³**), best-of-7 wall time, warmed `@allocated`. Moment
families as above (`V`/centroid/interface = `Tuple{0}`; `A`; `W`; `B`).

| Dim | Moment | Old time | New time | Speedup | Old alloc | New alloc | Reduction |
|---|---|--:|--:|--:|--:|--:|--:|
| 2D | `V`/centroid/interface | 51.8 ms | **14.4 ms** | ≈ 3.6× | 143.8 MB | **1.81 MB** | ≈ 79× |
| 2D | `A` | 75.8 ms | **57.0 ms** | ≈ 1.3× | 9.88 MB | **2.84 MB** | ≈ 3.5× |
| 2D | `W` | 90.1 ms | **32.4 ms** | ≈ 2.8× | 292.5 MB | **0.51 MB** | ≈ **572×** |
| 2D | `B` | 67.5 ms | **58.4 ms** | ≈ 1.2× | 6.77 MB | **3.03 MB** | ≈ 2.2× |
| 3D | `V`/centroid/interface | 2793 ms | **2380 ms** | ≈ 1.2× | 181.7 MB | **4.26 MB** | ≈ 43× |
| 3D | `A` | 5462 ms | **1599 ms** | ≈ 3.4× | 3.83 GB | **15.7 MB** | ≈ **244×** |
| 3D | `W` | 13772 ms | **4651 ms** | ≈ 3.0× | 11.76 GB | **5.10 MB** | ≈ **2198×** |
| 3D | `B` | 7564 ms | **1659 ms** | ≈ 4.6× | 11.82 GB | **5.08 MB** | ≈ **2219×** |

The 3D second-kind moments are the headline: `W` and `B` dropped from **≈ 11.8 GiB
to ≈ 5 MB per call** (≈ 2200×) while running 3–5× faster.

### Front path vs. signed-distance path

The two executors compute *different* objects — the front path returns exact
moments of a *given* discretized marker mesh (polygon/polyhedron clipping per cut
cell), while the level-set path returns moments of a *smooth implicit* surface
(analytic φ + the C/Julia VOFI engine). When a signed distance is available, the
implicit path is one to two orders of magnitude cheaper, since it never clips a
polytope. Same geometry/grids, ADELIE front vs. the best signed-distance backend
(`:vofi` level-set):

| Dim | Moment | ADELIE front | `:vofi` level set | Implicit advantage |
|---|---|--:|--:|--:|
| 2D | `V`/centroid/interface | 14.4 ms / 1.81 MB | **0.73 ms / 0.95 MB** | ≈ 20× faster |
| 3D | `V`/centroid/interface | 2380 ms / 4.26 MB | **12.3 ms / 2.67 MB** | ≈ 194× faster |
| 3D | `A` | 1599 ms / 15.7 MB | **3.60 ms / 1.09 MB** | ≈ 444× faster |
| 3D | `W` | 4651 ms / 5.10 MB | **43.1 ms / 1.15 MB** | ≈ 108× faster |

The gap is the cost of polytope clipping vs. analytic evaluation — expected, and
the reason the two paths coexist: use `:vofi`/`:vofijul` when you have an SDF, the
front path when you only have markers.

## Runtime

Eliminating per-cell heap traffic removes the garbage-collector pressure that
dominated the original hot loops. For `VofiJul`, the allocation collapse
translated into measured per-cell speedups of roughly **25–50×** (largest in 2D,
smallest in 4D, where the recursive slice path keeps more residual work):

| Path | Speedup (allocation-driven) |
|---|--:|
| 2D volume fraction | ≈ 50× |
| 3D volume fraction | ≈ 33× |
| 4D hypervolume fraction | ≈ 25× |

At the suite level, the `VofiJul` 4D test set (the allocation-heaviest path)
dropped from **≈ 1 min 53 s to ≈ 35 s**.

## Multithreading

Every geometry backend is **reentrant**: scratch is per-call, per-task, or
per-thread (no mutable module-global state), so user code can wrap its own
`Threads.@threads` over cells safely. This is locked by a "concurrent == serial"
test in each package's suite that runs under `-t` in CI.

| Package | Threaded path | Determinism |
|---|---|---|
| `VofiJul` | `vofi_get_cc_batch`, `vofi_get_cell_type_batch` | order-independent → bit-for-bit == serial |
| `isoap` | `isoapgrid` (fans cut cells across threads) | indexed write-back → bit-for-bit == serial at any thread count |
| `VOFTools` | user-driven (`@threads` over cells) | reentrant; concurrent == serial verified |

Reductions (e.g. assembling a moment field by `sparse(I,J,V)`) are written
order-independently so results never depend on thread scheduling.

## Capability matrix

What each backend gained relative to its original (the original codes were
`Float64`-only, serial, CPU-only):

| Package | Type-stable | Alloc-free hot path | Threaded | AD (ForwardDiff) | GPU (KernelAbstractions) |
|---|:--:|:--:|:--:|:--:|:--:|
| `VofiJul` | ✅ | ✅ | ✅ | ✅ geometry + level-set adjoint | ✅ CPU backend (device port pending) |
| `FrontCartesianGeometry` | ✅ | ✅ | ✅ (`:static`) | ✅ front-vertex / shape derivative | 🟡 point-in-front classify only |
| `isoap` | ✅ | ✅ | ✅ | ✅ field + vertex + Float32 | ✅ device-ready (`isoap_cell!` stack-only) |
| `VOFTools` | ✅ | ✅ | ✅ (reentrant) | ✅ plane + vertex + closure-param | ⬜ not started |
| `Vofinit` | ✅ | ✅ | ✅ (reentrant) | ⛔ C-library binding | ⛔ C-library binding |

🟡 = partial. `FrontCartesianGeometry` ships a KernelAbstractions point-in-front
classification kernel (`is_inside_gpu`, winding-number in 2D / solid-angle in 3D,
CPU-validated), but the **moment integrator** (the cut-cell polytope clip) has no
device kernel — it is CPU-only (type-stable, alloc-free, threaded). Its full
`integrate` AD path is differentiable through the front vertices, gated on a
`<:Real`-element mesh from `FrontIntrinsicOps` (verified `d(V_cell)/dR` AD-vs-FD
≈ 4e-10 in 2D and 3D).

**AD** covers both *geometry/shape* sensitivity (feed `ForwardDiff.Dual`
coordinates — the whole pipeline, root-finder included, is element-type generic)
and *level-set-parameter* sensitivity (custom adjoint rules where applicable).
All AD paths are finite-difference-verified across 1D–4D.

!!! warning "GPU on-device status"
    The KernelAbstractions kernels run today on the `CPU()` backend (validated
    bit-for-bit). `isoap`'s `isoap_cell!` is genuinely stack-only and ready for a
    real CUDA/ROCm device. `VofiJul`'s kernel still builds a per-work-item
    workspace on the heap (≈ 20 KB), so its on-device path awaits an
    immutable-workspace rewrite. `VOFTools` has no GPU kernel yet.
