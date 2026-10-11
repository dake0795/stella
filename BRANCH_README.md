# Branch `feature/continuous-flow-shear`

Continuous-in-time (moving-wavevector) treatment of the perpendicular equilibrium
flow shear in stella, together with a few related fixes and options of the explicit
time advance. All new options are off by default; with the default input the code is
unchanged.

The branch is based on commit `2b8e269f2addd0baa5991057eafa022135e04498` of `master`.

## Background

With a perpendicular flow shear `gamma_E = g_exb * g_exbfac` the radial wavenumber of
a mode of fixed `ky` evolves as `kx(t) = kx(0) - gamma_E ky t`. stella's standard
(Hammett) treatment keeps the `kx` grid fixed and remaps each mode by one grid spacing
`dkx` whenever its accumulated shift exceeds half a grid spacing, i.e. every
`T_ky = |dkx / (gamma_E ky)|`, so that the wavenumber of a mode follows a staircase in
time. The algorithm added here instead treats `akx(ikx)` as a label and evaluates every
term of the gyrokinetic and field equations at the physical wavenumber

    kx_phys(iky, ikx) = akx(ikx) - gamma_E * aky(iky) * tau(iky),

where `tau` is the residual shear time of the current Runge-Kutta stage. After each
accepted step the committed `tau` (stored in `shift_state`) is advanced; rows whose
`tau` has passed `T_ky/2` are relabelled by one grid spacing, which is the same physical
mode and therefore leaves the solution unchanged (apart from the mode that leaves the
grid). `kperp2`, the Bessel functions and the field denominators are recalculated
whenever the physical wavenumbers change. The approach is the one used for GS2 by
Christen, Barnes and Parra, J. Plasma Phys. 87, 905870230 (2021), "Continuous-in-time
approach to flow shear in a linearly implicit local delta-f gyrokinetic code"; the relabelling and the treatment of the
implicit terms are specific to stella.

## Commits

1. `Initialise flow-shear shift times when the shear rate is zero` -- bug fix:
   `shift_times` was uninitialised for `prp_shear_enabled = T, g_exb = 0` and entered
   the CFL estimate.
2. `Synchronise twist-and-shift junction points in the explicit schemes` -- option
   `explicit_junction_sync` (`'none'`, `'step'`, `'stage'`); default `'none'` is the
   original behaviour. Independent of the moving wavevector; the duplicated junction
   points of linked chains otherwise evolve independently in fully explicit
   simulations.
3. `Allow the flow-shear remap in fully explicit simulations` -- option
   `flow_shear_explicit_remap`; the standard remap is only applied in the implicit
   time advance.
4. `Split implicit mirror and streaming from the explicit flow-shear remap` -- with
   the above option, the implicit mirror and streaming are applied as a Lie split after
   the explicit step.
5. `Add a moving-wavevector treatment of the perpendicular flow shear` -- the feature
   (`flow_shear_moving_k`), module `gk_flow_shear_moving_k`, linear and fully explicit.
6. `Support the ExB nonlinearity with the moving-wavevector flow shear`.
7. `Support the implicit mirror with the moving-wavevector flow shear`.
8. `Support implicit parallel streaming with the moving-wavevector flow shear` --
   defect correction of the response-matrix solve to the physical wavenumbers.
9. `Add sub-cell response matrix sets as streaming preconditioners`.
10. `Share the flow-shear remap clock with processors without restart files` -- bug
    fix: processors without (kx,ky,z) points kept `shift_state = 0` after a restart,
    although they remap their part of g in the (vpa,mu,s) layout.
11. `Save and check the moving-wavevector phase in restart files`.
12. `Do not advance the flow-shear remap clock without shear` -- bug fix: the standard
    algorithm advanced `shift_state` with `g_exb = 0`, which a later restart with shear
    interpreted as an accumulated shift.
13. `Create the netcdf output file of a restarted simulation if it is missing` -- a
    simulation started from the restart files of another one stopped when opening
    its own output file.
14. `Add an optional free-energy budget diagnostic for explicit simulations`.

## Namelist `&flow_shear`: new inputs and defaults

| Input | Default | Meaning |
|---|---|---|
| `explicit_junction_sync` | `'none'` | Re-synchronise the duplicated junction points of linked modes in the explicit schemes: `'step'` after each accepted step, `'stage'` also before each Runge-Kutta stage. Required (`'step'` or `'stage'`) by the moving wavevector. |
| `flow_shear_explicit_remap` | `.false.` | Apply the fixed-grid remap after each accepted step of a fully explicit simulation. |
| `flow_shear_moving_k` | `.false.` | Moving-wavevector treatment of the flow shear. |
| `flow_shear_moving_k_import_zero_phase_seed` | `.false.` | Allow a restart file written without `flow_shear_moving_k`, at zero remap phase, to be used as the initial condition of a moving-k simulation. |
| `moving_k_stream_defect_tol` | `1.0e-10` | Relative field residual at which the defect correction of the implicit streaming stops. |
| `moving_k_stream_defect_maxit` | `30` | Maximum number of defect-correction iterations (the simulation stops if exceeded). |
| `moving_k_stream_defect_method` | `'gmres'` | `'simple'` (label response as approximate inverse), `'gmres'` (restarted GMRES preconditioned with it) or `'auto'` (`'simple'`, switching to `'gmres'` after 8 iterations). |
| `moving_k_stream_offset_sets` | `0` | Number `M` (even) of additional streaming response matrices at the sub-cell phases `m T_ky / M`, used as preconditioners for the lowest `ky` rows. |
| `moving_k_stream_offset_rows` | `3` | Number of non-zonal `ky` rows with additional response sets. |
| `moving_k_stream_offset_keep_on_dt_change` | `.false.` | Keep the additional sets as (stale) preconditioners when the time step changes instead of recalculating them. |
| `energy_budget_diagnostics` | `.false.` | Write the free-energy budget of the explicit time advance to `<run_name>.budget`. |

Restrictions of `flow_shear_moving_k` (checked at initialisation): local flux tube,
box grid with `nakx >= 3`, kinetic electrons, `phi` and `apar` without `bpar`,
`explicit_algorithm = 'rk4'`, `flip_flop = F`, `hammett_flow_shear = T`,
`omprimfac = 0`, no collisions, no sources, `D_hyper = 0` (parallel hyper-dissipation
is allowed), `explicit_junction_sync /= 'none'`. The mirror and the parallel streaming
may be implicit (applied as a split after the explicit step); the drifts may not. The
ExB nonlinearity requires `mod(nx, 3) /= 1`.

A typical input:

    &flow_shear
      prp_shear_enabled = .true.
      g_exb = 0.6
      omprimfac = 0.0
      flow_shear_moving_k = .true.
      explicit_junction_sync = 'stage'
      moving_k_stream_defect_method = 'auto'
      moving_k_stream_defect_tol = 1e-8
      moving_k_stream_offset_sets = 8
      moving_k_stream_offset_rows = 21
    /
    &numerical_algorithms
      explicit_algorithm = 'rk4'
      flip_flop = .false.
      stream_implicit = .true.
      mirror_implicit = .true.
    /

Restart files written with `flow_shear_moving_k` contain the variables `moving_k`,
`moving_k_convention`, `akx_label`, `aky_label` and `shear_rate` in addition to
`shiftstate` (the committed phase); a restart checks that they agree with the input.

## Tests

`AUTOMATIC_TESTS/numerical_tests/test_9_flow_shear/` adds to the standard test suite:
time traces against expected output for a linear and a nonlinear moving-k simulation and
for the explicit remap; the bookkeeping of the committed phase in the restart file; a
restarted versus an uninterrupted simulation; the independence of the defect correction
from its method, tolerance and preconditioner; the equality of the moving-k and standard
algorithms at zero shear; and the convergence of the time stepping.

## Verification status

The series was extracted from a development build of the same algorithm. Checks made on
this branch (binary built with `make ... SPLIT_BY_NUMA=`, i.e. one shared-memory domain
per node, for both the base commit and the branch):
- default input: output bitwise identical to the base commit on small grids with 1 and
  2 processors, for linear and nonlinear, electrostatic and electromagnetic simulations,
  with and without the fixed-grid remap in the implicit time advance, and for fully
  explicit sheared simulations; the outputs of the existing automated tests are bitwise
  identical as well. At 256 processors (a linear electromagnetic simulation with the
  fixed-grid remap and implicit streaming and mirror) the branch and the base commit
  differ by 4e-11 relative in the fields after the first time step, and the difference
  does not grow. Each binary is reproducible, and the intermediate commits of the
  branch, most of which do not touch this code path, alternate between the two results
  in no particular order. The cause has not been identified; it is also present without
  optimisation (DEBUG builds);
- every moving-wavevector configuration tested (implicit mirror, implicit streaming with
  and without additional response sets, fully explicit, nonlinear, linear, zero shear,
  energy budget, explicit-remap twin) bitwise identical to the development build at 1 and
  2 processors, and for a production-resolution linear simulation at 256 processors;
- restart equal to the uninterrupted simulation (bitwise at 1 and 2 processors); at 256
  processors the development build showed restart differences at the 1e-11 relative level,
  which it also showed for the standard algorithm with implicit streaming;
- zero shear with `flow_shear_moving_k = T` equals the fixed-grid explicit remap with the
  same junction synchronisation (bitwise for the linear test, 2e-10 for a nonlinear
  deck because of the order of the Fourier transforms); it differs from the standard
  implicit time advance at the 1e-5 level because that one does not synchronise the
  junction points during the explicit stages;
- the defect-corrected implicit streaming converges to the same solution with the gmres
  and auto methods and with or without additional response sets; the 'simple' method
  stalls at a relative residual of about 1e-10 and should be used with a tolerance of
  1e-8 or looser;
- the free-energy budget closes to the time-discretisation remainder (relative to the
  change of free energy: 5e-6 for the moving wavevector, 3e-6 for the explicit remap, in
  a bounds-checked nonlinear run with dt = 0.002);
- a bounds-checked (DEBUG) build runs the moving-wavevector configurations without
  errors and adds no compiler warnings to those of the base commit;
- the fully explicit moving-wavevector time stepping converges with the time step (rk4;
  the implicit mirror and streaming splits are first order in time).

Established on the development build and not repeated here: agreement of the sheared-frame
nonlinear term with an independent convolution; method-of-manufactured-solutions
convergence orders of the moving-wavevector operators; linear Floquet decay rates
independent of the radial grid spacing.

Not validated:
- long-time nonlinear behaviour of sheared turbulence with the moving wavevector
  compared with the fixed-grid treatment is under investigation; the continuous
  treatment is known to give different saturated fluxes from the discrete remap for
  some parameters (see Christen et al. 2021). Do not use for physics results until this
  is settled.

## Licence

stella is distributed under the MIT licence (`LICENCE.md`); the changes on this branch
are contributed under the same licence.
