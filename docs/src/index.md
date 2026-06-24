# ADELIE

**Adaptive Discretization Ecosystem for Liquid–Interface Evolution**

ADELIE is a collection of open-source Julia packages for the numerical simulation of multiphase flows, phase change, and moving interfaces. Packages are independently usable but designed to compose around a shared set of conventions.

## What we build

| Pillar | Focus |
|--------|-------|
| **Adaptive** | Structured and adaptive numerical methods that track interfaces without over-resolving the bulk |
| **Discretization** | Finite volume, cut-cell, and embedded-boundary schemes with clean component boundaries |
| **Ecosystem** | Modular, composable packages with a common design language |
| **Liquid–Interface** | Incompressible two-phase flows with surface tension, phase change, and topology changes |
| **Evolution** | Level-set, VOF, and hybrid interface-tracking methods for moving boundaries |

## Getting started

Browse the [Packages](packages.md) page to find the right starting point for your problem.

All packages are written in Julia and registered in the General registry. Install any package with:

```julia
using Pkg
Pkg.add("PackageName")
```

## Design principles

- Packages are independently usable but designed to compose
- Implementations target Julia for performance, expressiveness, and reproducibility
- Each package ships tests, documentation, and a versioned API following a shared standard

---

*ADELIE is named after the Adélie penguin — an animal that lives at the sharpest interface on Earth.*
