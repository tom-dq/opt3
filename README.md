# opt3

MayorFEM is a clean-slate nonlinear finite element prototype for macOS with a Metal element kernel and CPU fallback.

## Current Capabilities

- Nonlinear geometry via finite deformation gradient (`F`) in each tetrahedron.
- Nonlinear material with J2 (von Mises) return mapping and isotropic hardening.
- Simple ductility damage variable that softens stress as equivalent plastic strain accumulates.
- Enforced displacements (Dirichlet boundary conditions) with displacement-controlled load stepping.
- Newton solve with sparse global assembly and BiCGSTAB iterative linear solves.
- Element-level consistent tangent linearization (local directional derivatives through constitutive update).
- Metal backend for element residual/state evaluation; CPU backend remains available as fallback.
- Literature-style regression benchmark suite for yield transition and translation invariance checks.
- Stepwise load ramp snapshots for gradual loading analysis.
- Visualization export to VTK (`.vtk` + `series.pvd`) including stress, exaggerated deformation, and enforced displacement fields.

## Run

```bash
swift run mayor-fem --steps 12 --disp 0.08 --backend auto
```

Arguments:

- `--steps N` number of load steps
- `--disp value` end displacement on loaded face (x direction)
- `--backend auto|metal|cpu`
- `--steps N` controls gradual load increase from 0 to full prescribed displacement

### Visualization Export

```bash
swift run mayor-fem --steps 12 --disp 0.08 \
  --visualize out/vtk --deformation-scale 10 --backend auto
```

This writes one VTK file per converged load step plus `series.pvd` for ParaView time-series loading.

Included fields:

- Cell data: `von_mises`, `eq_plastic_strain`, `damage`, `load_factor`
- Point data: `displacement`, `prescribed_displacement`, `enforced_dof_count`, `prescribed_node`

Run benchmark suite:

```bash
swift run mayor-fem --benchmarks --backend cpu
```

## Tests

```bash
swift test
```

## Literature Anchors

This implementation follows standard building blocks from nonlinear solid mechanics and computational plasticity:

- Belytschko, Liu, Moran, and Elkhodary, *Nonlinear Finite Elements for Continua and Structures*.
- Simo and Hughes, *Computational Inelasticity*.
- Bonet and Wood, *Nonlinear Continuum Mechanics for Finite Element Analysis*.

The current code is tuned for small-model reliability and speed on laptop-scale runs; next steps are larger benchmark coverage and more element/integration variants.
