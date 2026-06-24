# CLAUDE.md

## Overview

**ADELIE-org.github.io** is the org-level landing site for the ADELIE organisation, built with [Documenter.jl](https://documenter.juliadocs.org). It has no Julia package of its own — it is a pure documentation site that catalogues all ADELIE packages.

## Repository layout

```
docs/
  make.jl          — Documenter build script
  Project.toml     — deps: Documenter only
  Manifest.toml    — generated, do not edit by hand
  src/
    index.md       — org overview and design principles
    packages.md    — one section per ADELIE package
.github/workflows/
  Docs.yml         — builds and deploys to gh-pages on push to main
```

## Commands

```bash
# Instantiate the docs environment (first time or after Manifest deletion)
julia --project=docs -e 'using Pkg; Pkg.instantiate()'

# Build the site locally
julia --project=docs docs/make.jl

# Output lives in docs/build/ — open docs/build/index.html in a browser to preview
```

## Adding a new package

1. Add a section to `docs/src/packages.md` following the existing pattern:
   - H3 heading linking to the GitHub repo
   - Badges (docs, CI)
   - One-paragraph description
   - Minimal code example
   - Install line
2. If the package belongs to a new category, add an H2 section heading.
3. Rebuild locally to verify no broken links or build warnings.

## Deployment

The `Docs.yml` workflow deploys to the `gh-pages` branch on every push to `main`. The `DOCUMENTER_KEY` secret must be set on the repo (generate with `DocumenterTools.genkeys()`). GitHub Pages must be configured to serve from `gh-pages`.

## Conventions

- No Julia package at the repo root — `docs/Project.toml` lists only `Documenter`.
- `inventory_version = "0"` is set in `make.jl` to suppress the version-extraction warning.
- `warnonly = true` keeps the build from failing on cross-reference issues during early development.
- Do not add a root `Project.toml` or `src/` — this is not a package.
