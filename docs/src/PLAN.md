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
   backends), `advance!`, velocity extension, and the full motion taxonomy of §4: four
   **regimes** (fixed, prescribed, rigid body, free boundary) and, for free boundary,
   the representation × updater cube — `ExplicitAdvection` wired to the classical
   advectors (**LevelSetMethods.jl** for SDF, **gVOF.jl** for VOF,
   **FrontTrackingMethods.jl** for front) and `InterfaceNewton` on an interface-balance
   residual with representation-specific corrections (**GlobalHeightFunctions** for SDF,
   **FrontCutTracking** Levenberg–Marquardt markers for front, VOF-field later).
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

### Path × capability — what each path can and cannot do

The coupling (§3), motion (§4) and nonlinearity (§5) machinery is **not uniformly
available on all three paths**, because they rest on two *different* Jacobian notions:

- **Solution-Jacobian via `tangent`** — ∂F/∂u, the operator differentiated w.r.t. its own
  field argument, as a reusable operator/matrix. This is the Newton and monolithic-
  cross-block enabler. **`:assemble` (assembled J) and `:matrixfree` (J·v action) only —
  `:kron` has no `tangent`.**
- **Parameter-AD** — ForwardDiff Duals flowing through moments → derivatives w.r.t.
  *parameters* (interface position, coefficients). **All three paths flow Duals** (verified),
  so shape/interface sensitivity and SciMLSensitivity work everywhere.

| Capability | `:kron` | `:assemble` | `:matrixfree` |
|---|:--:|:--:|:--:|
| Linear operators (G/H/div/lap/conv) | ✓ | ✓ | ✓ |
| Assembled matrix for direct solve (KLU/LU) | ✓ | ✓ | ✗ (matvec only) |
| Krylov / **Jacobian-free** Newton (J·v) | matvec | matvec | ✓ (native, GPU) |
| **Newton** solution-Jacobian (`tangent` ∂F/∂u) | ✗ | ✓ | ✓ |
| **Monolithic** coupled cross-blocks (§3) | ✗ | ✓ | ✓ |
| **Picard** / lagged linearization | ✓ | ✓ | ✓ |
| Parameter-AD (shape/interface sens., §4) | ✓ | ✓ | ✓ |
| Alloc-free time-loop refill (`update!`) | ✗ (rebuild) | ✓ | re-substitute |
| Role | oracle | workhorse | Krylov/GPU/large |

So: **`:kron` stays the linear-operator validation oracle** — it does Picard and
parameter-AD, but *not* `tangent`-based Newton or monolithic coupling, and it cannot
refill alloc-free in a time loop. All the nonlinear/monolithic/Newton machinery of
§3/§4/§5 runs on **`:assemble`** (direct solves, moving-geometry time loops) and
**`:matrixfree`** (Krylov, GPU, reverse-AD, large systems). `OneWay`/`Partitioned`
coupling is path-agnostic (each block picks its own path); only `Monolithic` and Newton
require assemble/matrixfree.

## 3. Coupling — field↔field multiphysics axes

Coupling is its own design surface, orthogonal to the physics packages and to time
stepping. A `CoupledProblem` is a set of **sub-problems** (each = operators +
conditions + unknowns produced by a physics package) plus a **`CouplingScheme`** that
says how their residuals/solves are combined. It lives in AdelieSolverCore; every
physics package emits sub-problems that plug in unchanged — no physics package knows
how it is coupled. Three orthogonal axes:

### Axis 1 — coupling topology (how strongly the blocks talk)

- **One-way** (`OneWay`): sub-problem B reads A's converged solution as data; A never
  sees B. No iteration, single ordered pass. E.g. flow → passive scalar transport.
- **Two-way / partitioned** (`Partitioned`): each block keeps its own system; iterate
  to consistency within the step/steady solve. Sub-strategies: block Gauss–Seidel
  (sequential, latest values) / block Jacobi (parallel, lagged); fixed-point (Picard)
  with optional Aitken / IQN relaxation for stiff coupling; or Newton on the coupling
  residual — **Jacobian from FieldAlgebra `tangent`, analytic, not FD**.
- **Monolithic** (`Monolithic`): one `LinearSystem`/`NonlinearSystem` over all
  unknowns; off-diagonal coupling blocks assembled through the shared layout; single
  solve. The cross-blocks are `tangent`s of one physics' residual w.r.t. the other's
  unknowns — free from the symbolic layer. **This is the thing Penguin could not do
  cheaply** (it hand-coded coupling blocks per pair).

Accuracy ladder: one-way < partitioned-lagged < partitioned-iterated-to-tol ≈
monolithic; cost ladder is the reverse. Partitioned-iterated and monolithic reach the
*same* fixed point — monolithic gets there in fewer (superlinear) iterations at the
cost of assembling/solving the full block system.

### Axis 2 — solution context (where the coupling loop lives)

- **Steady** (`SteadyCoupling`): no time loop; the coupling iteration *is* the outer
  solver — fixed-point or Newton over the coupled residual to `‖r‖ < tol`. One-way
  degenerates to a single ordered pass.
- **Unsteady** (`TimeLoop`): coupling is resolved every step, inside the θ/SciML
  stepper. This is the only context where splitting (Axis 3) applies.

### Axis 3 — time-splitting (unsteady only: how operators share a step)

- **Unsplit** (`Unsplit`): all physics advanced together each step (monolithic-in-time,
  or a within-step partitioned iteration). No splitting error.
- **Lie / Godunov** (`Lie`): sequential sub-steps A then B over Δt. First-order in Δt,
  cheapest.
- **Strang** (`Strang`): ½Δt A, Δt B, ½Δt A. Second-order, symmetric — the default
  when a split is wanted.
- **IMEX / additive** (`IMEX`): implicit on the stiff term (diffusion, pressure),
  explicit on the non-stiff (advection, reaction) — a tableau, not a sequential split.
  For advection–diffusion–reaction with a stiff diffusion or reaction.
- **Iterated split** (`IteratedStrang` / predictor–corrector): repeat the split
  sub-steps and iterate to tol → recovers the unsplit fixed point, trading splitting
  error for iterations. The bridge between partitioned and monolithic *within* a step.

Splitting error and coupling-iteration error are **distinct**: an unsplit-but-partitioned
step still carries partitioned-iteration error until converged; a Strang-split-but-
monolithic-per-substep step still carries O(Δt²) splitting error. Convergence studies
report both.

### How the axes compose (and the interface-motion special case)

A run picks one value per axis: `(Monolithic, TimeLoop, Unsplit)` for tightly-coupled
Boussinesq; `(OneWay, TimeLoop, Lie)` for passive transport; `(Partitioned+Newton,
SteadyCoupling, —)` for steady FSI.

**Decided: field↔field and field↔geometry coupling stay separate abstractions.**
The interface-motion coupling of Phase 5 (`ExplicitAdvection` / `InterfaceNewton`,
geometry↔field) is kept as its own machinery, *not* folded into the `CouplingScheme`
above. The two share vocabulary and read as analogues —

- `ExplicitAdvection` is the geometry↔field analogue of `(OneWay-lagged, Lie)`;
- `InterfaceNewton` is the analogue of `Monolithic` / `Partitioned`-Newton —

but the interface residual (geometry moves, moments rebuild, DOF layout changes shape)
is different enough from a fixed-layout field block that forcing one abstraction would
distort both. `CouplingScheme` governs **field↔field** coupling only; motion coupling
lives in InterfaceMotion (Phase 5). Revisit unification only if the two converge in
practice.

### API sketch

```julia
prob = CoupledProblem(
    (:flow => flow_sub, :heat => heat_sub),               # named sub-problems
    couplings = [(:flow, :heat) => BoussinesqCoupling()], # who reads whom + the cross-term
    scheme    = Monolithic(),                             # OneWay() | Partitioned(Picard; relax=Aitken) | Monolithic()
)
solve!(prob, SteadyCoupling(; tol = 1e-8))               # steady: coupling iteration = outer solve
step!(prob, TimeLoop(; θ = 0.5, split = Strang()), Δt)   # unsteady: coupling resolved per step
```

Axis 1 (`scheme`) and Axes 2/3 (time driver + `split`) are independent arguments;
`CoupledProblem` validates legal combinations (`split` is illegal under
`SteadyCoupling`; `Monolithic` scheme + `Strang` split = monolithic solve *within* each
split sub-step).

## 4. Interface motion — the geometry↔field coupling

The other coupling surface: how the *interface itself* moves. Kept as its own machinery,
**separate from the §3 field↔field `CouplingScheme`** (decided). Organized by **motion
regime**, in increasing order of coupling strength. The regime — not a scheme flag —
decides what machinery is active.

### Regime 1 — Fixed
Geometry static. No slab; `SpaceMoments` only. The baseline for Phases 1–4.

### Regime 2 — Prescribed
Interface position is known a priori as a function of time. The representation carries a
**time component directly** — `SpaceTimeLevelSet` (SDF with time), `STFront`, `STPLIC` —
so `integrate_spacetime` produces `SpaceTimeMoments` (slab: Vn, Vnp1, tγ) and the
space-time operators + solver do the rest. **No feedback, no updater, no coupling loop —
it is automatic**: prescribing motion = fixing the interface-velocity input to the
space-time predictor (Penguin's `solve_prescribed_motion!` reused the exact slab
machinery with an externally supplied `Vprescribed`). The cheap, fully-determined regime.

### Regime 3 — Rigid body
Interface moves as a **rigid transform** driven by a **6-DOF ODE** (Newton–Euler) on the
fluid forces (+ gravity, contact/collision). The coupling is **ODE↔PDE**: fluid solve
↔ body-state ODE.
- The body *shape* is rigid, but the geometry still moves each step, so **space-time
  slab moments must be recomputed every step** — the swept slab between tⁿ and tⁿ⁺¹
  (Vn, Vnp1, tγ) depends on the motion, via the same `integrate_spacetime` path as
  prescribed motion. What is cheap is producing the interface *position* (a rigid
  transform of the body — rotate/translate the markers or shift the SDF argument, no
  shape re-reconstruction); the slab moment integration over the moving geometry is
  **not** skipped.
- Staggered-explicit first (predict body state over [t,t+dt] → set moving geometry →
  fluid solve → integrate hydro load: force + torque via momentum/embedded-force balance
  → advance the 6-DOF ODE; Penguin `NavierFSIProblem`/`MultiBodyFSIProblem`, symplectic
  Euler). Strong (Picard/Newton) coupling later for stiff FSI (light/neutrally-buoyant
  bodies).
- **Contact/collision is new** (Penguin had none): enters as a force/constraint in the ODE.

Spiritually this is a §3-style sub-problem (a body-DOF block) more than a free-boundary
interface residual — but it lives here under motion, driven by the geometry layer.

### Regime 4 — Free boundary
Interface velocity comes **from the fields** (Stefan flux, kinematic condition, GCL, FSI
traction). This is where the **representation × updater** cube lives. Port PenguinStefan's
`AbstractInterfaceRep` (it is good): `phi_values`, `predict_phi`, `advance!`,
`extend_speed!`, `coupling_mode`.

**Interface representation** (3):
- **Level-set / SDF** (`LevelSetRep`) — φ on the grid; ADELIE geometry via VofiJul/implicit.
- **Front tracking** (`FrontTrackingRep`) — Lagrangian markers (`CurveMesh`/`SurfaceMesh`/
  `PointFront1D`, FrontIntrinsicOps).
- **VOF** — volume-fraction field + PLIC.

**Interface updater** (2 strategies, geometry↔field coupling):
- **`ExplicitAdvection`** — classical lagged advection of the representation by the
  extended interface speed. One classical-advection backend per representation:
  - SDF → **LevelSetMethods.jl** (`LevelSetEquation`, `AdvectionTerm`, `integrate!`,
    `reinitialize!`, `extend_along_normals!`).
  - VOF → **gVOF.jl** (`faceflux!` + `vofadv!`, `reconstruct!` PLIC).
  - Front → **FrontTrackingMethods.jl** (`FrontEquation`, `integrate!`,
    `redistribute!`/`remesh!`).
  Cheap; CFL-limited accuracy in the coupling.
- **`InterfaceNewton`** — inner correction: Newton iterations per step on an **interface
  balance residual** (thermal Stefan, species Stefan, GCL, Gibbs–Thomson, kinematic).
  The residual drives a **representation-specific update** of the interface DOFs:
  - SDF → update a **global height function** (GlobalHeightFunctions: `update_heights`,
    residual `Rcol = (Hⁿ⁺¹−Hⁿ) + Fcol/ρL`), rebuild φ via `phi_from_xf` — Penguin
    `:ghf_newton` / `:ls_newton`.
  - Front → update **markers via Levenberg–Marquardt** on marker normal displacements
    (FrontCutTracking intercept Jacobian `J[k,I]=ρL·A_{k,I}`; LM `A = JᵀWJ + R + λD`,
    `Δd = A\(−g)`, curvature/mass regularization, backtracking on λ) — Penguin `:ft_lm`.
  - VOF → update the **VOF field** (roadmap; Penguin had no VOF-Newton).
  **ADELIE advantage:** the residual Jacobian comes **analytically via parameter-AD
  through the slab moments** (verified ~1e-10 vs FD for level-set and STFront), replacing
  Penguin's finite-difference volume Jacobians. STPLIC is discontinuous in α
  (plane-dilation classification) → VOF defaults to `ExplicitAdvection` (Newton only with
  a smoothed classification, if ever needed).
- **Explicit seeds Newton:** `ExplicitAdvection` provides the predictor / initial guess
  for `InterfaceNewton` each step (Penguin does exactly this).

Not every cube cell is needed; the **regime** selects the machinery: Fixed/Prescribed
need no updater; Rigid body is the 6-DOF ODE coupling; Free boundary is the only regime
where representation × updater matters.

## 5. Nonlinearity — one residual/Jacobian, Newton everywhere

Nonlinearity is cross-cutting, not a per-physics afterthought. Sources across the target
physics:

- **Nonlinear advection** — NS convection u·∇u; advection whose velocity depends on the
  transported field.
- **Reaction** — Arrhenius R(T, Yₖ), rate laws; often stiff.
- **Solution-dependent coefficients** — κ(T), μ(T), ρ(T,Yₖ); harmonic/arithmetic
  averaging of a solution-dependent diffusivity.
- **Nonlinear conditions** — Robin / radiation (T⁴), Gibbs–Thomson (curvature), the
  Stefan flux balance.
- **EOS / constraint** — low-Mach ρ = P₀/RT and the heat-release divergence constraint.
- **Geometry** — free-boundary moments depend nonlinearly on interface position (§4).

**The core enabler.** Every ADELIE operator is a symbolic FieldAlgebra tree, and
`tangent` returns its **exact Jacobian-as-stencil analytically** — as a matrix-free
`apply`/`pullback` (path 1) or an assembled CSC (path 2). So a nonlinear residual `F(x)`
built from the same operators yields `J = tangent(F)` for free: **no hand-coded Newton
blocks (Penguin's pain), no finite-difference Jacobians, no black-box AD of an opaque
assembly.** This is the *same* primitive that produces the coupling cross-blocks (§3) and
the interface-residual Jacobian (§4). One abstraction, three uses.

**Solve strategies** (in AdelieSolverCore, over `LinearSolve` / `NonlinearSolve`):

- **Picard / lagged linearization** — freeze the nonlinear factor (velocity → Oseen,
  coefficient → previous iterate), solve linear, repeat. Robust, linear rate; the cheap
  default and the standard NS workhorse. Anderson acceleration optional.
- **Newton** — `F(x)=0`, `J = tangent(F)`, solve `J δ = −F`, update. Quadratic; the
  default when Picard stalls or high accuracy is needed. Cross-term Jacobians (reaction
  ∂R/∂Yⱼ, convection ∂(u·∇φ)/∂u) come from `tangent`.
- **Jacobian-free Newton–Krylov** — use matrix-free `apply(J,·)` (path 1) inside a Krylov
  solve for large/coupled systems; never assemble J (`apply` == exact J·v, verified).
- **Modified / quasi-Newton** — freeze and reuse the J factorization across iterations/steps
  (the `LinearSystem` cached factorization already supports this); refactor only when
  convergence degrades.
- **Globalization** — line search / damping, and continuation/homotopy (ramp Re, ramp
  reaction rate, ramp Gibbs–Thomson) for stiff or far-from-solution starts. The
  Picard/explicit result seeds Newton (mirrors §4's explicit-seeds-Newton).
- **Stiffness** — IMEX (§3 Axis 3) treats the stiff nonlinear reaction implicitly (its own
  Newton on the reaction block) while advection stays explicit; or go fully implicit
  Monolithic-Newton.

**Where nonlinearity meets coupling — the nesting.** The *same* Newton/Picard machinery
runs at three scopes, and the §3 topology axis decides how they nest:

- **intra-block** — nonlinearity inside one physics (NS convection, reaction);
- **inter-block** — the coupled multiphysics residual (§3);
- **geometry** — the interface residual (§4).

`Monolithic` fuses all active scopes into **one Newton** over the full unknown vector
(one `J = tangent(F_total)`, superlinear, tightest). `Partitioned` keeps them as **nested
loops** — inner Newton per block, outer Picard/Newton over the coupling, outermost
interface correction. Same residual/Jacobian primitive at every level; the choice is cost
vs. robustness, exactly as in §3. A run picks, per scope, Picard-vs-Newton and, across
scopes, monolithic-vs-nested — independently.

## 6. Time integration — steppers, slabs, and the mass matrix

Time stepping is a mechanism in its own right, and it is where coupling (§3),
nonlinearity (§5), and the moving-geometry moments meet. Two paradigms, deliberately
kept distinct by **whether the geometry moves within the step**:

### Paradigm A — Method of lines (MOL), for fixed geometry
Discretize space (the §1–§2 operators) → a semi-discrete ODE `M u̇ = −A(u) u + b(u,t)`
→ hand to a time integrator. The **mass matrix `M` is the cut-cell capacity `V`**
(diagonal); the spatial operator is constant in time, so this is clean and lets us
inherit SciML's adaptive, stiff, and IMEX integrators. Interface-DOF elimination reduces
the DAE to an `ODEProblem` with mass matrix (Phase 6 MOL export); steady problems become a
`NonlinearProblem`.

### Paradigm B — Space-time / slab, for moving geometry
When the geometry moves *within* the step (prescribed / rigid / free-boundary, §4), the
spatial operator is **not** constant over `[tⁿ, tⁿ⁺¹]`, so MOL does not hold. Instead
integrate the **space-time slab** directly, using `SpaceTimeMoments` (Vn, Vnp1, tγ): the
step is a two-level slab, not a line. This is why ADELIE carries space-time moments at
all. The **geometric conservation law (GCL)** is satisfied *by construction* — the
discrete volume change `Vnp1 − Vn` matches the swept-interface flux carried by `tγ`, so a
uniform field is preserved under mesh motion (no spurious source). Prescribed motion is
the pure form (slab known a priori); rigid/free-boundary rebuild the slab each step (§4).

### The stepper hierarchy (what lands when)
1. **θ-method** one-step, hand-rolled: explicit (θ=0), Crank–Nicolson (θ=½, default),
   implicit Euler (θ=1). Reuses the cached `LinearSystem` factorization when dt and
   operators are unchanged. The workhorse for Phases 2–5, both paradigms.
2. **SciML integrators** (Phase 6, fixed-geometry MOL): SDIRK/ESDIRK, Rosenbrock(-W),
   BDF/FBDF for stiff, explicit RK for non-stiff — with error-controlled **adaptive dt**.
   IMEX (`SplitODEProblem`; KenCarp/ARS tableaus) for the stiff-reaction / non-stiff-
   transport split of §3 Axis 3.
3. **Slab stepper** (moving geometry): the θ or a 2-level slab scheme over the
   space-time moments; dt is CFL- or Newton-controlled (no SciML adaptivity here yet).

### How the step composes with Newton (§5) and coupling (§3)
An implicit step *defines the per-step residual* `F(uⁿ⁺¹) = 0`; §5 solves it (Picard or
Newton), and §3 + Axis-3 splitting organize the multiphysics *inside* the step:
- **Fully implicit, unsplit, monolithic** = one Newton per step over all unknowns
  (`J = tangent(F_total)`, §5) — tightest, on `:assemble`/`:matrixfree`.
- **IMEX / split** = implicit (stiff) terms enter the step's Newton; explicit terms are
  evaluated at `uⁿ`. `Strang`/`Lie` sequence the sub-steps (§3 Axis 3).
- **Partitioned** = each block does its own implicit solve; an outer Picard/Newton over
  the coupling closes the step.
The stepper is thus orthogonal to §3/§5: it fixes *what* residual each step poses; those
sections fix *how* it is solved and nested.

### Path dependence
θ-stepper + direct solve rides `:assemble` (cached factorization, alloc-free `update!`
refill in the loop — the moving-geometry hot path); Krylov/JFNK stepping rides
`:matrixfree`; `:kron` is not used in production time loops (no alloc-free refill — §2).

**Open items:** MOL-vs-slab is chosen by geometry (fixed→MOL/SciML, moving→slab), not by
the user; adaptive dt for moving geometry (error control through a changing mass matrix)
is deferred. See Open decision on SciML depth (§10).

## 7. Physics coverage — mapping the goals onto the layers

The target physics — **N scalars per phase, flow+thermal, low-Mach combustion, and
(later) triphasic** — are combinations of a small set of primitives, *not* new
infrastructure. Two knobs carry almost all of it: the **layout** (AdelieSolverCore's
`UnknownLayout`, parametric in phase count × scalars per phase) and the **§3 coupling
schemes**. What each goal needs:

- **N scalars per phase** — the base case of AdelieScalar. `UnknownLayout(phases,
  scalars)`; each scalar has its own coefficients / source / border+interface conditions;
  inter-scalar coupling (reactions) is a block over scalars with the reaction Jacobian
  from `tangent` (`Monolithic`) or a `Partitioned` sweep. No new package. (Phase 3.)

- **Flow + thermal** (forced / natural convection, Boussinesq, Rayleigh–Bénard) — a
  `CoupledProblem` of AdelieFlow (u,p) + AdelieScalar (T): `OneWay` for forced convection
  (flow→T), `Monolithic` / `Partitioned`-Newton when buoyancy feeds back (T→flow via the
  Boussinesq body force). This is the §3 API-sketch example verbatim — no infrastructure
  beyond the coupling term. (Phase 7, on top of Phases 4 + 3.)

- **Low-Mach combustion** — the demanding goal; the only one needing genuinely new
  *model* pieces (still no new infrastructure). It couples flow ↔ thermal ↔ N species ↔
  stiff chemistry with **variable density**:
  - AdelieFlow gains a **variable-density, constrained-divergence** mode: `∇·u = Q` set
    by heat release / species change (not zero), thermodynamic pressure P₀(t) split from
    the dynamic pressure, EOS `ρ = P₀/(RT)` closing the loop.
  - AdelieScalar carries **N species + temperature** with a **stiff Arrhenius reaction
    source** → `IMEX` (implicit chemistry/diffusion, explicit transport) or a `Monolithic`
    Newton on the reaction block (`tangent` gives the chemistry Jacobian — §5).
  - The whole thing is one `CoupledProblem` (flow + thermal + species) resolved
    `Monolithic` or tight-Picard per step, with `Strang`/`IMEX` splitting to isolate the
    stiff chemistry. So low-Mach = {AdelieFlow variable-density extension} + {AdelieScalar
    N-species+reaction} + {§3 coupling with IMEX}. Schedule after flow and scalar are both
    solid (Phase 8 territory).

- **Triphasic (later)** — falls out of the **parametric layout**: `UnknownLayout(3, …)`
  adds a third phase block plus the extra interface/jump rows (N-phase conditions from
  CutCellConditions). The math and solver ride the same code exercised for mono/diph;
  only the **geometry is genuinely harder** — triple junctions / three-region moments
  need the geometry layer to return per-phase moments at junction cells (a reconstruction
  problem, flagged, not a solver problem). Generality lives in the type now; the
  triphasic *implementation* lands when a triphasic geometry backend does.

Throughline: **every goal = layout(phases × scalars) + a §3 coupling recipe + (for
low-Mach) a variable-density flow model.** Block assembly, AD Jacobians, and the solver
are shared — the modularity Penguin lacked.

## 8. Phases

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
   Picard/Newton drivers, `CoupledProblem` (see §3). Phase 2 delivers the coupling
   skeleton: Axis 1 `OneWay` + `Partitioned(Picard)` + a `Monolithic` slot; Axis 2
   `SteadyCoupling` + `TimeLoop` drivers; Axis 3 `Unsplit` + `Lie`. Enough to run a
   two-scalar reaction–diffusion coupled case end-to-end.
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
- Coupling (§3) matured for scalars: Axis 3 `Strang` and `IMEX` split drivers
  (stiff diffusion/reaction implicit, advection explicit); inter-scalar reaction
  coupling available `Monolithic` via `tangent` or `Partitioned` per scalar block.
- Moving geometry for all of the above (prescribed motion first).
- Gate: parity vs PenguinDiffusion/PenguinTransportDiffusion + convergence suites.

### Phase 4 — flow
- Staggered cut-cell layouts (velocity components on face grids, pressure on cells) —
  moments per component grid already come from the geometry layer.
- Darcy first (smallest saddle-point), then Stokes (gauges: pin/mean pressure),
  Brinkman, then Navier–Stokes (Picard + Newton via symbolic `tangent` of the
  convection term — replaces Penguin's hand-coded Newton blocks; see §5).
- Gate: Penguin parity on Stokes lid-driven/Poiseuille + NS benchmark (drag/lift
  on cylinder), diphasic Stokes jump case.

### Phase 5 — motion and phase change

The interface-motion taxonomy (regimes Fixed / Prescribed / Rigid body / Free boundary,
and the free-boundary representation × updater cube) is specified in **§4**. This phase
implements it, in InterfaceMotion.jl + AdeliePhaseChange.jl, roughly in regime order:

- **InterfaceMotion.jl**: the `AbstractInterfaceRep` abstraction (level-set / front / VOF
  backends) with `phi_values`/`predict_phi`/`advance!`/`extend_speed!`; the
  `ExplicitAdvection` updater wired to LevelSetMethods.jl (SDF), gVOF.jl (VOF),
  FrontTrackingMethods.jl (front); the `InterfaceNewton` updater with the height-function
  (SDF), Levenberg–Marquardt marker (front) and VOF-field (roadmap) corrections, all on
  **AD residual Jacobians** through the slab moments. Prescribed motion needs no updater
  (space-time SDF → moments, §4 Regime 2); rigid-body 6-DOF ODE coupling (§4 Regime 3),
  staggered-explicit first.
- **AdeliePhaseChange.jl**: Stefan mono → diph, stepper = predict geometry → slab moments →
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
- Coupling maturity (§3): complete the `Monolithic` scheme (multi-physics blocks
  assembled into one system via layouts); add `Partitioned`-Newton with Aitken/IQN
  relaxation and `IteratedStrang`; and SciML splitting (`SplitODEProblem` / IMEX
  tableaus) alongside the MOL export.
- FSI/collision, low-Mach combustion, triphasic, elasticity ride on this foundation
  (each is a new model + coupling recipe, not new infrastructure — see §7 for how the
  physics goals map onto layout + coupling, and §5 for the nonlinear solve; low-Mach also
  needs the variable-density flow model).

## 9. Cross-cutting rules

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

## 10. Open decisions (flagged, with recommendation)

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
5. **Splitting driver** (§3 Axis 3): hand-rolled θ + split drivers (recommended for
   Phases 2–3 — small, full control, AD-transparent) vs SciML `SplitODEProblem` /
   OrdinaryDiffEq IMEX from the start. Recommendation: hand-rolled first, adopt SciML
   splitting in Phase 6 with the MOL export, once the interfaces are stable.
6. **Coupling default**: symbolic `tangent`-based cross-blocks make `Monolithic` and
   `Partitioned`-Newton cheap, so they are first-class from Phase 2 (not deferred as in
   Penguin). `OneWay`/`Partitioned(Picard)` remain the low-cost defaults; `Monolithic`
   is opt-in per `CoupledProblem` where the coupling is stiff.
7. **Nonlinear driver** (§5): NonlinearSolve.jl for Newton with the **analytic `tangent`
   Jacobian supplied** (recommended — not FD, not AD-of-black-box), plus a hand-rolled
   Picard/Anderson loop for the cheap lagged cases; vs a fully hand-rolled Newton. Reuse
   the `LinearSystem` cached factorization for modified-Newton. Recommendation: Picard
   available from Phase 2 (needed for NS anyway), NonlinearSolve Newton wired in when the
   first genuinely nonlinear physics lands (advection–diffusion reaction / Phase 3, flow /
   Phase 4).
