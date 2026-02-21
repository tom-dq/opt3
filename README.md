# opt3

MayorFEM is a 2D nonlinear finite element prototype focused on small-model speed and reliability for laptop workflows.

## Current Capabilities

- 2D nonlinear geometry (finite deformation kinematics in plane strain embedding).
- Nonlinear material with J2 (von Mises) plasticity, isotropic hardening, and ductile damage softening.
- Enforced displacements with gradual load stepping (displacement-control ramp).
- Explicit nonlinear solver path tuned for small-model robustness.
- CPU and Metal backends for the explicit element loop (`auto` selects Metal when available).
- Implicit Newton path retained for comparison/debugging.
- Triangular element workflow:
  - Linear triangles (`tri3`)
  - Quadratic triangles (`tri6`) represented and solved through consistent sub-triangulation
- Uniform mesh subdivision (`--subdivide`) for refinement studies.
- Patch-based topology optimization:
  - Density design variables per element
  - Reference element toggled on/off in local patches
  - Arbitrary patch objective callback API
  - Compliance objective provided out of the box
- Visualization export per load step:
  - VTK series (`.vtk` + `series.pvd`)
  - PNG frame series
  - Includes stress, exaggerated deformation, and enforced displacement overlays.

## Run (Explicit)

```bash
swift run mayor-fem --solver explicit --steps 12 --disp 0.08 --nx 12 --ny 4 --order linear
```

Key arguments:

- `--solver explicit|implicit`
- `--steps N` number of load increments
- `--disp value` final prescribed displacement on loaded edge
- `--backend auto|cpu|metal`
- `--nx N`, `--ny N` base rectangular mesh divisions
- `--order linear|quadratic` element order
- `--subdivide L` uniform refinement levels (each level multiplies triangle count by 4)
- Explicit controls:
  - `--explicit-substeps N`
  - `--relax-iters N`
  - `--dt value`
  - `--damping value`
  - `--mass-density value`
  - `--density-penalty value`

### Visualization Export

```bash
swift run mayor-fem \
  --solver explicit --steps 12 --disp 0.08 --order quadratic --subdivide 1 \
  --visualize out/viz2d --deformation-scale 10 \
  --image-width 1600 --image-height 1000
```

This writes:

- `out/viz2d/vtk/step_XXXX.vtk` + `out/viz2d/vtk/series.pvd`
- `out/viz2d/png/step_XXXX.png`
- `out/viz2d/densities.csv`

Included fields:

- Cell data: `von_mises`, `eq_plastic_strain`, `damage`, `load_factor`, `density`
- Point data: `displacement`, `prescribed_displacement`, `enforced_dof_count`, `prescribed_node`

### Benchmarks (Literature-Anchored)

```bash
swift run mayor-fem --benchmarks --solver explicit --backend cpu
```

Current benchmark checks include:

- Below-yield elastic response (Simo/Hughes-style anchor)
- Above-yield plastic response (Belytschko-style loading path)
- Subdivision consistency
- Translational invariance (Bonet/Wood-style patch anchor)

### Topology Optimization

```bash
swift run mayor-fem \
  --topopt --solver explicit --backend metal \
  --steps 8 --disp 0.05 --nx 6 --ny 2 \
  --topopt-iters 8 --patch-radius 1 \
  --min-density 0.05 --max-density 1.0 \
  --move-limit 0.12 --reference-stride 1 \
  --objective compliance \
  --visualize out/topopt
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
