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
| `isoap` | ✅ | ✅ | ✅ | ✅ field + vertex + Float32 | ✅ device-ready (`isoap_cell!` stack-only) |
| `VOFTools` | ✅ | ✅ | ✅ (reentrant) | ✅ plane + vertex + closure-param | ⬜ not started |
| `Vofinit` | ✅ | ✅ | ✅ (reentrant) | ⛔ C-library binding | ⛔ C-library binding |

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
