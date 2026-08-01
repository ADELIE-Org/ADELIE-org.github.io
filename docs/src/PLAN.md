# ADELIE implementation plan — from geometry to multiphysics

**Status: revised 2026-08-01.** Supersedes the 2026-07-02 plan, which was written
when the stack stopped at the algebra layer. Phases 0–5 of that plan have landed:
operators, conditions, solver core, the scalar physics package, interface motion,
flow, Darcy, and monophasic phase change are all implemented, registered, and
gated. This revision records **what exists**, then plans the remaining work
against the three target goals — **multicomponent phase change, combustion, FSI**.

Sections §3–§6 of the old plan (coupling axes, motion regimes, nonlinearity,
time integration) described a design that is now built. They are kept here,
re-cast in the present tense, because the design is still the contract.

## 0. Where we stand

### The stack, bottom to top

| Layer | Packages | Status |
|---|---|---|
| Grids | `CartesianGrids` | ✅ `CartesianGrid{N,T}`, `SpaceTimeCartesianGrid{N,T}` |
| Geometry | `CartesianGeometry` (router), `VofiJul`, `Vofinit`/`VOFTools`, `FrontCartesianGeometry`, `FrontIntrinsicOps`, `FrontSpaceTimeDEC`, `isoap` | ✅ `SpaceMoments`/`SpaceTimeMoments` for **all three** interface reps (level-set, PLIC/VOF, front), static + space-time slab, curvature, AD/GPU/threaded |
| Algebra | `AlgebraCore`, `ScalarAlgebra`, `FieldAlgebra` | ✅ symbolic fields, `δ₊/δ₋`, `tangent`, matrix-free `apply`/`pullback`, `assemble`/`assemble!` |
| Operators | `AdelieOperators` | ✅ G/H/div/grad/laplacian/convection on three paths; `MomentBundle` adapter; moving operators; mobility; staggered flow blocks |
| Conditions | `CutCellConditions` | ✅ borders (Dirichlet/Neumann/Robin/Periodic), mono interface closures, diphasic jump rows, velocity/traction jumps (incl. strain form), Darcy rows |
| Solver core | `AdelieSolverCore` | ✅ `UnknownLayout`, `BlockSystem`, `SystemOperator`, `newton!`/`picard!`, `ThetaIntegrator`, `MovingIntegrator` + mono/diph/Darcy moving slabs, `CoupledProblem` (OneWay/Partitioned/Monolithic), `CoupledTimeLoop`, SciML DAE bridge |
| Physics | `AdelieScalar`, `AdelieFlow`, `AdelieDarcy`, `AdeliePhaseChange`, `InterfaceMotion` | partial — see below |
| Validation | `AdelieValidation` | ✅ case groups A–G with convergence gates |

### Physics coverage, honestly

- **`AdelieScalar`** — effectively complete. Mono/diphasic × N species × {steady,
  transient, prescribed-moving} × {diffusion, advection–diffusion, reaction
  (Newton, analytic or AD Jacobian), `κ(x)`, `κ(u)` (Picard)}, plus `dae_problem`
  for every family and time-aware sources throughout. **This is the
  multicomponent-species infrastructure** the combustion and phase-change goals
  need; it is already built.
- **`AdelieFlow`** (v0.8) — Stokes and Navier–Stokes on the staggered cut-cell
  MAC grid, steady and unsteady, monolithic (no projection, no splitting);
  diphasic Stokes/NS with a **fixed** interface, velocity continuity + prescribed
  normal-traction jump or `σκ` surface tension; `interface_force`/
  `interface_torque`; Boussinesq coupling both partitioned and monolithic.
  Validated: canonical Stokes (inclined channel, Basilisk couette/wannier,
  Womersley), Schäfer–Turek 2D-1/2D-2, Hadamard–Rybczynski, de Vahl Davis to
  Ra = 1e6. **Fixed geometry only** — see §1.
- **`AdelieDarcy`** — mixed saddle Darcy, mono/diphasic, tensor mobility,
  gravity, interface force, and **prescribed moving geometry** via the slab.
- **`AdeliePhaseChange`** — **monophasic** Stefan only. Explicit (level-set and
  front), `:ghf_newton` (digit-for-digit Penguin parity), `:ft_lm`
  (~38–83× the explicit accuracy at equal `Δt`). Frank disk/sphere, 1-D Neumann
  similarity, planar MMS all gated.
- **`InterfaceMotion`** — `LevelSetRep`, `GlobalHFRep`, `FrontTrackingRep`;
  `ExplicitAdvection` + `InterfaceNewton` (`:ghf_newton`, `:ft_lm`). No VOF
  representation yet.

## 1. The gap to the three goals

### Multicomponent phase change

Monophasic Stefan works. Missing:

- **diphasic Stefan** — two-phase temperature across a moving interface. Core
  already has `moving_diph_slab`, so this is mostly `AdeliePhaseChange` work.
- **Gibbs–Thomson** — the curvature it needs is already produced by the geometry
  layer (`CartesianGeometry.curvature`); the condition type is unwritten.
- **species-driven interface motion** — a mass-flux Stefan condition, equilibrium
  / Henry partitioning at a *moving* interface, composition-dependent `Tm`.
- **a joint `(T, Y⃗)` interface residual** for `InterfaceNewton` — today the
  residual is scalar (thermal only).
- **Stefan flow** — the volume expansion at the interface needs moving-geometry
  flow (below).

### Combustion

**Nothing exists.** A repo-wide search for `low-mach | combustion | arrhenius |
variable.density | reacting` returns no implementation hits. What exists is the
reusable half: N-species transport with a Newton-coupled reaction source and
analytic cross-Jacobians. Missing:

- **variable-density, constrained-divergence flow**: `∇·u = Q` set by heat
  release and species change, thermodynamic pressure `P₀(t)` split from the
  dynamic pressure, EOS `ρ = P₀/(RT)` closing the loop.
- **chemistry**: an Arrhenius reaction model, `μ(T)`, mixture-averaged species
  diffusion.
- **IMEX** — `AdelieSolverCore`'s `IMEX` splitting is *defined but throws*. Stiff
  chemistry needs it (or a fully monolithic Newton on the reaction block).

### FSI

The **load side is done**: `interface_force`/`interface_torque` are
Hadamard-checked, with action–reaction agreement to 1e-12. That is the hard half
of the fluid→body transfer. Missing:

- **moving-geometry flow** (below) — the actual blocker.
- **the 6-DOF Newton–Euler body**, rigid-transform interface motion, and the
  body-state↔fluid coupling loop (§4 Regime 3). Staggered-explicit first, then
  strong coupling for light / neutrally-buoyant bodies where added mass makes the
  explicit scheme unstable.
- **contact / collision** — new; Penguin had none.
- **elastic FSI is out of reach** — it needs a solid-mechanics package that does
  not exist. That is a scope decision, not a task.

### The single blocker under all three

**`AdelieFlow` has no moving-geometry path.** `AdelieScalar` and `AdelieDarcy`
both ride `integrate_spacetime` → space-time slab; `AdelieFlow` does not — there
is no `moving`/`slab` code in its `src/`. FSI needs it directly; Stefan flow and
boiling need it; any moving-interface combustion needs it. It is also the place
where discrete GCL and degenerate small-cell rows bite hardest, so it should be
budgeted as genuinely hard, not as a port.

## 2. Architecture (built — this is the contract)

### The three operator paths

An operator is *defined once*, symbolically, in FieldAlgebra; the path is a
materialization choice, and all three live behind one API in `AdelieOperators`.

| Capability | `:kron` | `:assemble` | `:matrixfree` |
|---|:--:|:--:|:--:|
| Linear operators (G/H/div/lap/conv) | ✓ | ✓ | ✓ |
| Assembled matrix for direct solve (KLU) | ✓ | ✓ | ✗ (matvec only) |
| Krylov / Jacobian-free Newton (J·v) | matvec | matvec | ✓ (native, GPU) |
| Newton solution-Jacobian (`tangent` ∂F/∂u) | ✗ | ✓ | ✓ |
| Monolithic coupled cross-blocks (§3) | ✗ | ✓ | ✓ |
| Picard / lagged linearization | ✓ | ✓ | ✓ |
| Parameter-AD (shape/interface sensitivity) | ✓ | ✓ | ✓ |
| Alloc-free time-loop refill (`update!`) | ✗ (rebuild) | ✓ | re-substitute |
| Role | oracle | workhorse | Krylov/GPU/large |

`:kron` is the **linear-operator validation oracle**, not a solve path — BC and
jump rows break the tensor structure anyway. All Newton / monolithic / moving
machinery runs on `:assemble` (direct, the moving hot path) and `:matrixfree`
(Krylov, GPU, large systems).

### Coupling — field↔field (`AdelieSolverCore`)

A `CoupledProblem` is a set of named `SubProblem`s plus a `CouplingScheme`. No
physics package knows how it is coupled. Three orthogonal axes:

1. **Topology** — `OneWay` (B reads A's converged solution) / `Partitioned`
   (block Gauss–Seidel fixed point, constant or Aitken relaxation) /
   `Monolithic` (one system over all unknowns; cross-blocks are `tangent`s of one
   physics' residual w.r.t. another's unknowns — the thing Penguin could not do
   cheaply). Accuracy ladder one-way < partitioned-lagged <
   partitioned-to-tol ≈ monolithic; cost ladder reversed.
2. **Solution context** — `SteadyCoupling` (the coupling iteration *is* the outer
   solve) / `TimeLoop` (resolved every step).
3. **Time splitting** (unsteady only) — `Unsplit` / `Lie` / `Strang` /
   `IteratedStrang` / `IMEX`. **`IMEX` throws today**; `IteratedStrang` is
   currently ≡ `Unsplit`.

Splitting error and coupling-iteration error are distinct; `step_report!` reports
the latter numerically, the former still needs a reference solve.

**Decided and unchanged: field↔field and field↔geometry coupling stay separate
abstractions.** `CouplingScheme` governs field↔field only; interface motion lives
in `InterfaceMotion` (§4). They read as analogues — `ExplicitAdvection` ≈
`(OneWay-lagged, Lie)`, `InterfaceNewton` ≈ `Monolithic` — but the interface
residual moves the geometry and reshapes the DOF layout, which is different
enough that one abstraction would distort both.

### Interface motion — field↔geometry (`InterfaceMotion`)

Organized by **regime**; the regime, not a flag, selects the machinery.

1. **Fixed** — static geometry, `SpaceMoments` only.
2. **Prescribed** — the representation carries time directly
   (`SpaceTimeLevelSet`/`STFront`/`STPLIC`), `integrate_spacetime` produces the
   slab, no updater and no coupling loop. ✅ implemented for scalar and Darcy.
3. **Rigid body** — a 6-DOF Newton–Euler ODE driven by the fluid loads. The body
   *shape* is rigid but the geometry still moves, so slab moments are rebuilt
   every step; what is cheap is producing the interface *position*. ❌ **not
   implemented** (see §3 P2).
4. **Free boundary** — the interface velocity comes from the fields. This is
   where the **representation × updater** cube lives:

   | | `ExplicitAdvection` | `InterfaceNewton` |
   |---|---|---|
   | Level-set (`LevelSetRep`) | ✅ LevelSetMethods.jl | ✅ `:ghf_newton` (global height function) |
   | Front (`FrontTrackingRep`) | ✅ marker displacement | ✅ `:ft_lm` (Levenberg–Marquardt on normal displacements) |
   | VOF | ❌ planned (gVOF.jl) | ✗ by design — STPLIC is discontinuous in `α` |

   Explicit always seeds Newton as the per-step predictor.

### Nonlinearity — one residual/Jacobian primitive

Every ADELIE operator is a symbolic FieldAlgebra tree and `tangent` returns its
exact Jacobian-as-stencil — assembled or matrix-free. So a nonlinear residual
built from the same operators yields its Jacobian for free: no hand-coded Newton
blocks, no finite-difference Jacobians, no black-box AD of an opaque assembly.
The *same* primitive gives the coupling cross-blocks and the interface-residual
Jacobian. One abstraction, three uses.

The machinery runs at three scopes — **intra-block** (NS convection, reaction),
**inter-block** (coupled multiphysics), **geometry** (the interface residual) —
and the topology axis decides the nesting: `Monolithic` fuses the active scopes
into one Newton; `Partitioned` keeps them as nested loops. `newton!`/`picard!`
reuse a cached factorization (`reuse=true`) on the fixed sparsity pattern.

Still missing here: **line search / trust region** in `newton!` (full steps
only), and JFNK preconditioning.

### Time integration — two paradigms, chosen by the geometry

- **Paradigm A — method of lines**, for **fixed** geometry: the mass matrix is
  the cut-cell capacity `V` (diagonal, singular on γ rows ⇒ an index-1 DAE), the
  spatial operator is constant, so SciML's adaptive/stiff/IMEX integrators apply.
  `dae_problem` covers every scalar family and flow.
- **Paradigm B — space-time slab**, for **moving** geometry: the operator is not
  constant over `[tⁿ, tⁿ⁺¹]`, so MOL does not hold. The step integrates the slab
  using `SpaceTimeMoments`. **GCL holds by construction** — the discrete volume
  change `Vⁿ⁺¹ − Vⁿ` matches the swept-interface flux carried by `tγ`. Ψ⁺/Ψ⁻
  multiply **on the right** of K/C (flux of the time-blended field), which is
  what keeps flux antisymmetry at fresh/dead cells.

Adaptive `Δt` for moving geometry goes through SciML as a genuine implicit DAE
(`moving_dae_problem`), not a hand-rolled controller.

## 3. Roadmap

Each item ends with a hard validation gate and the 5-objective checklist
(type-stable / alloc / threads / GPU / AD).

### P1 — Moving-geometry flow (the keystone)

`moving_stokes_slab` / `moving_ns_slab` in `AdelieSolverCore`, mirroring the
existing `moving_mono_slab`/`moving_diph_slab`, plus the `AdelieFlow` driver and
a `MovingLevelSet` interface spec (as `AdelieScalar` and `AdelieDarcy` already
have).

**Gates:** free-stream preservation on a translating cut wall (the discrete GCL
test — this is the one that fails first); prescribed oscillating cylinder against
a body-fitted reference; then moving *diphasic* with `σκ` — a rising-bubble case.

**Budget it as hard.** The two things that decide whether the order survives are
the space-time capacities and the linearised traction row; degenerate
small-volume fragment rows are the failure mode to watch, not the small cells
themselves.

### P2 — FSI (rides P1 immediately)

`InterfaceMotion` Regime 3: 6-DOF Newton–Euler state + rigid transform of the
representation. Staggered-explicit coupling first (predict body state → set
moving geometry → fluid solve → integrate the hydrodynamic load → advance the
ODE), then `Partitioned`-Newton with Aitken (already in Core) for light bodies
where added mass destabilizes the explicit scheme. Contact enters as a force /
constraint in the ODE.

**Gates:** settling / rising cylinder terminal velocity; vortex-induced
vibration; a rigid Turek–Hron variant. Elastic FSI is explicitly **out of scope**
until a solid-mechanics package is decided on.

### P3 — Diphasic Stefan + Gibbs–Thomson

Mostly `AdeliePhaseChange` work — Core's `moving_diph_slab` is already there.
Two-phase interface flux balance; a `GibbsThomson` condition in
`CutCellConditions` fed by `CartesianGeometry.curvature`; per-phase `κ`.

**Gates:** two-phase Frank sphere; a Gibbs–Thomson dendrite-tip case; explicit
vs `:ghf_newton` vs `:ft_lm` at equal `Δt`.

### P4 — Multicomponent phase change

The moving diphasic N-scalar machinery already exists in `AdelieScalar`. The new
pieces are the species interface condition (equilibrium / Henry partitioning
**with** mass transfer), the joint `(T, Y⃗)` interface residual for
`InterfaceNewton`, and composition-dependent `Tm`. Stefan flow lands here via P1.

**Gates:** binary-alloy solidification with solute rejection; an evaporating
droplet with a species-driven `ṁ`.

### P5 — Low-Mach combustion

- `AdelieFlow` gains the **variable-density constrained-divergence** mode:
  `∇·u = Q`, `P₀(t)` split from the dynamic pressure, EOS `ρ = P₀/(RT)`.
- `AdelieScalar` gains **Arrhenius chemistry** on the existing N-species reaction
  block (`tangent` supplies the chemistry Jacobian).
- `AdelieSolverCore` **implements `IMEX`** (currently throws) — this needs
  sub-problems that expose an explicit stiff/non-stiff *term* split, a richer
  abstraction than today's `TransientSubProblem`.
- The whole thing is one `CoupledProblem` (flow + thermal + species) resolved
  `Monolithic` or tight-Picard per step. **The monolithic Boussinesq validation
  is the exact structural template** for the flow⇄T⇄species Newton.

**Gates:** 1-D premixed flame speed vs a reference mechanism; a counterflow
diffusion flame.

### Ordering note

**P1+P2 and P3+P4 are independent.** P3 and most of P4 need no moving flow —
pure-conduction multicomponent Stefan is reachable without touching
`AdelieFlow`. If the near-term objective is a multicomponent-phase-change result,
do **P3+P4 first** and let P1 wait.

## 4. Cross-cutting debts

- **`AdelieSolverCore`**: `IMEX` throws; `IteratedStrang` ≡ `Unsplit`;
  `Partitioned`-Newton missing (Picard only); `Monolithic` requires the user to
  hand-fuse sub-problems into one `BlockSystem` — a layout-stitching helper would
  remove that manual step; no line search in `newton!`; no threading; no GPU path
  (`SystemOperator` is GPU-ready *in principle*, untested, no smoke test).
- **Triphasic** is blocked upstream, not in the solver: `UnknownLayout` /
  `BlockSystem` already stack arbitrary named blocks, but `jump_closure` ties
  exactly **two** phases. Triple-junction jump rows in `CutCellConditions` come
  first, and the geometry layer must return per-phase moments at junction cells —
  a reconstruction problem, flagged, not a solver problem.
- **Scalar border convection**: `AdelieScalar` has no border convection closure,
  so a flow with `u·n ≠ 0` through the box borders carries an `O(h)` border
  truncation (measured 4e-3 at n=41). Recorded, not fixed.
- **VOF representation** in `InterfaceMotion` (explicit only — STPLIC is
  discontinuous in `α`).
- **Docs**: this site's `packages.md` still stops at the geometry layer — the
  entire operator/solver/physics stack (`AdelieOperators` … `AdeliePhaseChange`)
  is uncatalogued. `PLAN.md` is not in `make.jl`'s `pages` and `pagesonly = true`,
  so it is not built into the site at all.

## 5. Cross-cutting rules (unchanged)

- **Validation style**: element-by-element parity vs old Penguin wherever Penguin
  has the capability (this caught real bugs in the geometry and phase-change
  ports), plus manufactured solutions where it doesn't.
- **5 objectives** per package: JET/`@inferred` type-stability, `@allocated`
  gates on hot loops, `:static` threading, KernelAbstractions device paths,
  ForwardDiff Dual flow (both geometry-parameter and solution AD).
- **Linear algebra**: LinearSolve.jl exclusively; KLU default.
- **No physics package touches raw moments** — physics sees only
  `AdelieOperators` operators and `CutCellConditions` closures. Each physics
  package concentrates its moment access in a single `geometry.jl`.
- **Deps resolve through `ADELIERegistry`**, never `[sources] {path=...}`.

## 6. Open decisions

1. **Elastic FSI** — needs a solid-mechanics package (new) or an external
   coupling. Not scoped. Rigid-body FSI (P2) does not depend on this.
2. **Triphasic** — the layout is already parametric; the implementation waits on
   a triple-junction geometry backend *and* N-phase jump rows. Generality stays
   in the type, not in speculative implementation.
3. **Chemistry source** — hand-rolled Arrhenius vs a mechanism reader
   (Cantera-style). Recommendation: hand-rolled one-step and two-step first, so
   P5's flow and coupling work is not gated on a parser.
4. **GPU reach** — the geometry and operator layers are GPU-capable; the solver
   front door is KLU-direct. A GPU Krylov path would ride Core's matrix-free
   `SystemOperator`. Deferred until a problem size demands it.
