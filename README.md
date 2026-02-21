# opt3

MayorFEM is a clean-slate nonlinear finite element prototype for macOS with a Metal element kernel and CPU fallback.

## Current Capabilities

- Nonlinear geometry via finite deformation gradient (`F`) in each tetrahedron.
- Nonlinear material with J2 (von Mises) return mapping and isotropic hardening.
- Simple ductility damage variable that softens stress as equivalent plastic strain accumulates.
- Enforced displacements (Dirichlet boundary conditions) with displacement-controlled load stepping.
- Newton solve with finite-difference Jacobian for robust small-model convergence.
- Metal backend for element residual/state evaluation; CPU backend remains available as fallback.

## Run

```bash
swift run mayor-fem --steps 12 --disp 0.08 --backend auto
```

Arguments:

- `--steps N` number of load steps
- `--disp value` end displacement on loaded face (x direction)
- `--backend auto|metal|cpu`

## Tests

```bash
swift test
```

## Literature Anchors

This implementation follows standard building blocks from nonlinear solid mechanics and computational plasticity:

- Belytschko, Liu, Moran, and Elkhodary, *Nonlinear Finite Elements for Continua and Structures*.
- Simo and Hughes, *Computational Inelasticity*.
- Bonet and Wood, *Nonlinear Continuum Mechanics for Finite Element Analysis*.

The current code is a pragmatic prototype for small models; next steps for production-grade use are consistent tangents, sparse global assembly/solvers, and broader element/integration support.
