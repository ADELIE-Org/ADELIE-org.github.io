# ADELIE implementation plan — from geometry to multiphysics

Status of this document: proposed plan (2026-07-02). Grids + geometry + algebra are done;
this plan covers everything above them.

## 0. Where we stand

**Done (registered, 5-objective-hardened):**

| Layer | Packages | What they give the next layer |
|---|---|---|
| Grids | CartesianGrids | `CartesianGrid{N,T}`, `SpaceTimeCartesianGrid{N,T}`, node coords, spacing |
| Geometry | CartesianGeometry (router), VofiJul, Vofinit/VOFTools, FrontCartesianGeometry, FrontIntrinsicOps, FrontSpaceTimeDEC, isoap | `SpaceMoments{N,T}` / `SpaceTimeMoments{N,T}`: V, centroid, interface, cell_type, A, W, B (+ Vn/Vnp1/tγ for slabs), for **all three interface reps** (level-set, PLIC/VOF, front), static + space-time, AD/GPU/threaded |
| Algebra | AlgebraCore, ScalarAlgebra, FieldAlgebra | symbolic fields, `δ₊/δ₋`, `tangent` (Jacobian-as-stencil), matrix-free `apply`/`pullback`, `assemble`→SparseMatrixCSC, zero-alloc `assemble!` |

**The gap:** nothing yet connects moments to algebra. There are no named cut-cell
operators (G, H, div, grad), no boundary/interface conditions, no problem/solver layer,
no physics. Old Penguin Gen 2 (`CartesianOperators` → `PenguinSolverCore` → per-physics
models) is the architectural template to port — cleaned, deduplicated, and rebuilt on
FieldAlgebra.

**Discretization being ported (unchanged math):** capacities as diagonal weights,
`G[d] = δ₋∘B[d]`, `H[d] = A[d]·δ₋ − δ₋∘B[d]`, `Winv`, bulk/interface (ω/γ) block systems:
`[Iᴰ GᵀW꜀G  IᴰGᵀW꜀H ; Iᵦ HᵀW꜀G  Iᵦ HᵀW꜀H + Iₐ Γ]`, diphasic = two blocks + jump rows.

## 1. Target package map (new packages, dependency order)

Names are proposals; ADELIE precedent is to reuse Penguin-era names when the package is
the successor (as done for CartesianGeometry).

1. **CartesianOperators.jl** — the moment→operator bridge. The single most important
   next package.
2. **CutCellConditions.jl** (successor of PenguinBCs) — border BCs + interface/jump
   conditions, algebra-aware.
3. **AdelieSolverCore.jl** (successor of PenguinSolverCore) — layouts, LinearSystem,
   DOF masking, LinearSolve/NonlinearSolve integration, θ/SciML time stepping, coupling.
4. **AdelieScalar.jl** — diffusion, advection–diffusion, reaction–diffusion; mono/di/N-phase;
   1..N scalars; fixed + moving geometry.
5. **InterfaceMotion.jl** — the `AbstractInterfaceRep` abstraction (level-set, front, VOF
   backends), `advance!`, velocity extension, and the motion taxonomy: three **motion
   regimes** (prescribed, rigid body, free boundary) × two **interface updaters**
   (classical explicit advection; Newton inner-correction on an interface residual).
   (Absorbs PenguinVOFKinematics, LevelSetMethods glue, front advection, VCCIS/GHF later.)
6. **AdelieFlow.jl** — Darcy, Stokes, Navier–Stokes, Brinkman (staggered cut-cell
   saddle-point systems, gauges, Picard/Newton).
7. **AdeliePhaseChange.jl** — Stefan (mono/diph), Gibbs–Thomson, Stefan–Stokes,
   Darcy phase change, Boussinesq/Rayleigh–Bénard.
8. **AdelieIO.jl / AdelieViz.jl / AdelieAnalysis.jl** — VTK output, plotting,
   convergence tooling (ports of VTKOutputs, PenguinPlots, PenguinAnalysis). Low
   priority; port on demand.

Fewer, fatter physics packages than Penguin (which had ~15 solver repos with the block
assembly reimplemented in each). One scalar package, one flow package, one phase-change
package; variants are types/dispatch inside them, not new repos.

## 2. The three operator paths — how they coexist

All three paths live in CartesianOperators.jl behind one API. An operator is *defined
once*, symbolically, in FieldAlgebra; the path is an execution/materialization choice:

- **Path 1 — matrix-free:** `tangent` the symbolic form once → `FieldOuterSparse`
  stencil; `apply(J, x)` for J·x, `pullback(J, w)` for Jᵀ·w. For Krylov solves,
  GPU, and reverse-mode AD.
- **Path 2 — sparse CSC:** `assemble(J, ranges)` → `SparseMatrixCSC` for direct
  factorizations (LU/KLU); `assemble!` refills values in the time loop (already
  verified zero-alloc). This is the workhorse path, equivalent to what Penguin built.
- **Path 3 — Kronecker (legacy/reference):** direct port of Penguin
  `CartesianOperators` (`kron`-lifted `ẟ_m/δ_p/Σ_m/Σ_p` × diagonal capacities →
  explicit G/H/Winv sparse matrices). Kept primarily as the *validation oracle*: paths
  1/2 must reproduce its matrices element-by-element, exactly as CartesianGeometry was
  validated element-by-element against old Penguin. Once parity is locked in CI it can
  fade to a test dependency.

API sketch:

```julia
ops = CutCellOperators(moments::SpaceMoments, grid; path=:assemble)  # :matrixfree | :assemble | :kron
# named operators built from δ₊/δ₋ + bound moments:
gradient(ops), divergence(ops), interface_gradient(ops)   # G, -(Gᵀ+Hᵀ), H
laplacian(ops; coeff, mode=:harmonic)                     # GᵀW꜀G etc.
convection(ops, u; scheme=:upwind1)                       # C, K  (:centered/:upwind1/:advective/:skew)
update!(ops, moments)                                     # assemble!-based refill, alloc-free
```

## 3. Phases

Each phase ends with a hard validation gate and a 5-objective checklist
(type-stable / alloc / threads / GPU / AD). Nothing advances past a gate that isn't green.

### Phase 0 — hygiene (days)
- Fix FieldAlgebra docs drift: WIP.md claims `assemble`/`pullback` don't exist; they do
  and are tested. Reconcile WIP.md/AGENTS.md/README with the code before building on it.
- Write the moment↔algebra contract doc (one page): node-shaped AoS arrays bind via
  `substitute` as coefficient arrays; halo convention; `cell_type` masking; the 8 known
  reconciliation points from the geometry audit (cell-centred α vs node-shaped moments,
  flat executor outputs, backend-dependent `interface_centroid`, `SpaceMoments` vs
  `SpaceTimeMoments` membership, GPU-reduced output set, implicit fresh/dead masks).

### Phase 1 — CartesianOperators.jl (the next real step)
1. **Port the Kronecker path** from Penguin `CartesianOperators.jl` (types.jl,
   differential_operators.jl, assembly.jl, capacity.jl), consuming ADELIE
   `SpaceMoments` instead of Penguin capacities. This is quick and gives the oracle.
2. **Express the same operators in FieldAlgebra**: G/H/Winv/div/grad/laplacian/convection
   as symbolic trees over `δ₊/δ₋` with moment arrays bound by `substitute`; SVector
   moments (A, W, B) via `FieldRef` components. `tangent` → assemble (path 2) and
   apply/pullback (path 1).
3. **MomentBundle adapter**: one struct that normalizes `SpaceMoments` /
   `SpaceTimeMoments` / GPU-reduced outputs for the operator layer; derives fresh/dead
   masks (`Vn==0 & Vnp1>0` etc.); the space-time operator convention uses **Vn, Vnp1, tγ**
   (not Wt/Bt), matching the validated slab plan.
4. **Validation gate**: path 2 CSC == path 3 Kronecker matrices, element-by-element,
   1D/2D/3D, all three interface reps, static + one space-time slab case; path 1
   `apply` == mat-vec of the assembled matrix to machine precision.
5. **5-objective pass** (the algebra layer currently has none of threads/GPU):
   `assemble!` alloc-free in a rebuild loop with Dual-eltype moments (parameter AD
   end-to-end: interface param → moments → operator entries); threaded `assemble`
   band loops; KernelAbstractions kernel for matrix-free `apply` on device arrays.
   GPU can land as a follow-up objective within the phase, not a blocker for Phase 2.

### Phase 2 — conditions + solver core + first physics (the vertical slice)
1. **CutCellConditions.jl**: port PenguinBCs types (Dirichlet/Neumann/Robin/Periodic,
   ScalarJump/FluxJump/RobinJump/GibbsThomson, Inflow/Outflow/Traction/PressureOutlet…).
   Two application mechanisms, both supported: (a) symbolic — conditions contribute
   `Iₐ/Iᵦ`-style operator blocks and RHS terms before assembly (AD-transparent, preferred);
   (b) row stamping on assembled matrices (fallback, needed for exact Penguin parity).
2. **AdelieSolverCore.jl**: `UnknownLayout` (mono/diph → N-phase generalization),
   `LinearSystem{T}` with cached factorization, inactive-DOF masking (the
   `remove_zero_rows_cols` logic, done once, mask-based not copy-based),
   `solve!` via LinearSolve.jl only (drop IterativeSolvers), θ-stepper,
   Picard/Newton drivers, `CoupledProblem` (one-way / Picard two-way / monolithic slot).
3. **AdelieScalar.jl v0: monophasic diffusion**, steady + unsteady, variable
   coefficient (harmonic/arithmetic), on all three operator paths.
4. **Validation gate**: element-by-element parity of A and b vs Penguin
   `DiffusionUnsteadyMono` on identical geometry; manufactured-solution convergence
   (2nd order bulk); a moving-geometry prescribed-motion case via `SpaceTimeMoments`
   reproducing a Penguin prescribed-motion result.
5. This phase proves the whole architecture; expect design iteration here, budget for it.

### Phase 3 — scalar completion
- Diphasic diffusion (4×4 block + ScalarJump/FluxJump rows). Layout is parametric in
  (phase count, scalar count); mono + diph exercised now, triphasic later rides on the
  same layout without a rewrite.
- Advection–diffusion (schemes: centered/upwind1/advective/skew-symmetric), reaction
  terms, N scalars per phase — thermal, species, … — each with its own coefficients,
  sources, and border/interface conditions (block layout over scalars; inter-scalar
  coupling via reaction Jacobians — FieldAlgebra `tangent` gives these for free).
- Moving geometry for all of the above (prescribed motion first).
- Gate: parity vs PenguinDiffusion/PenguinTransportDiffusion + convergence suites.

### Phase 4 — flow
- Staggered cut-cell layouts (velocity components on face grids, pressure on cells) —
  moments per component grid already come from the geometry layer.
- Darcy first (smallest saddle-point), then Stokes (gauges: pin/mean pressure),
  Brinkman, then Navier–Stokes (Picard + Newton via symbolic `tangent` of the
  convection term — replaces Penguin's hand-coded Newton blocks).
- Gate: Penguin parity on Stokes lid-driven/Poiseuille + NS benchmark (drag/lift
  on cylinder), diphasic Stokes jump case.

### Phase 5 — motion and phase change

Three orthogonal axes, kept separate in the design:

- **Motion regime** — *prescribed* (geometry known a priori; slab moments per step,
  no feedback), *rigid body* (interface moves as a rigid transform driven by a 6-DOF
  ODE on fluid forces/collisions; geometry update is cheap — transform, don't
  re-reconstruct), *free boundary* (interface velocity comes from the fields:
  Stefan flux, kinematic condition, FSI traction balance).
- **Interface updater** — two strategies behind one interface:
  - `ExplicitAdvection`: classical lagged coupling — extend the speed, advance the
    representation, rebuild slab moments, step the field solve. Cheap, CFL-limited
    accuracy in the coupling.
  - `InterfaceNewton`: inner correction — Newton iterations per time step on an
    **interface residual** (e.g. flux/volume mismatch between the interface position
    and the field solution), converging geometry and fields together. Successor of
    Penguin's `:ghf_newton`/`:ls_newton` coupling modes and the front-tracking
    Gauss-Newton volume-Jacobian — but where Penguin used finite-difference Jacobians,
    ADELIE gets the residual Jacobian **analytically via parameter-AD through the slab
    moments** (already verified ~1e-10 vs FD for level-set and STFront). STPLIC is
    known discontinuous in α through plane-dilation classification, so the VOF backend
    defaults to `ExplicitAdvection` (Newton only with a smoothed classification, if
    ever needed).
- **Interface representation** — level-set, front (FrontIntrinsicOps meshes), VOF;
  port PenguinStefan's `AbstractInterfaceRep` design (it is good): `phi_values`,
  `predict_phi`, `advance!`, `extend_speed!`.

Not every cell of the regime × updater × representation cube is needed: prescribed
motion needs no updater; rigid body pairs with explicit force coupling first (strong
coupling later for FSI); free boundary is where both updaters matter.

- AdeliePhaseChange.jl: Stefan mono → diph, stepper = predict geometry → slab moments →
  scalar solve → interface speed from flux → advance (explicit) or correct (Newton);
  Gibbs–Thomson; then Stefan–Stokes, Darcy phase change, Boussinesq.
- Gate: Frank-sphere / plane-front analytic Stefan benchmarks with **both updaters**
  (Newton should beat explicit on accuracy at equal dt); a rigid-body prescribed-force
  case; Penguin parity where results exist.

### Phase 6 — SciML + coupling maturity
- MOL export: interface-DOF elimination → `ODEProblem` with mass matrix (port the
  PenguinMOL reduction onto AdelieSolverCore); `NonlinearProblem` for steady/Newton;
  SciMLSensitivity-compatible parameter AD (the operator layer is already
  Dual-transparent, so this is mostly plumbing + tests).
- Monolithic coupling slot: assemble multi-physics blocks into one system via layouts.
- FSI/collision, low-Mach combustion, triphasic, elasticity ride on this foundation
  (each is a new model + coupling recipe, not new infrastructure).

## 4. Cross-cutting rules

- **Validation style**: element-by-element parity vs old Penguin wherever Penguin has
  the capability (this caught real bugs in the geometry port), plus manufactured
  solutions where it doesn't.
- **5 objectives**: every package repeats the geometry-layer discipline — JET/`@inferred`
  type-stability tests, `@allocated` gates on hot loops (`update!`/`assemble!`/`apply`),
  `:static` threading with maxthreadid-safe pools, KernelAbstractions device paths,
  ForwardDiff Dual flow tests (both geometry-parameter AD and solution AD).
- **Linear algebra**: LinearSolve.jl exclusively; KLU default for the block systems;
  no IterativeSolvers dependency.
- **No physics package touches raw moments**: physics sees only CartesianOperators
  operators + CutCellConditions. That's the modularity boundary Penguin never had.

## 5. Open decisions (flagged, with recommendation)

1. **Package names** above are placeholders — reuse Penguin names where a package is a
   direct successor (recommended, matches CartesianGeometry precedent) or all-new
   ADELIE names.
2. **BC mechanism**: symbolic-first with stamping fallback (recommended) vs
   stamping-only parity port first. Symbolic-first costs more in Phase 2 but is what
   makes end-to-end AD and GPU work later.
3. **Phase/scalar layout** — *decided*: scope is monophasic and diphasic now,
   triphasic later; N scalars per phase (thermal, species, …) from the start. So
   `UnknownLayout` is parametric in (phase count, scalars per phase), with only 1- and
   2-phase paths exercised until triphasic physics arrives — generality in the type,
   not speculative implementation.
4. **SciML depth in Phase 2**: LinearSolve only (recommended) vs pulling ODEProblem
   in immediately. The θ-stepper port is small and keeps Phase 2 focused.
