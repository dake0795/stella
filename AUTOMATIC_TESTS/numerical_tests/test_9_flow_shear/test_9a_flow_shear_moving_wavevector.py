################################################################################
#          Check the moving-wavevector treatment of the flow shear             #
################################################################################
# Perpendicular flow shear with a moving wavevector (flow_shear_moving_k) and
# the related options of the explicit time advance (explicit_junction_sync,
# flow_shear_explicit_remap). The flags we test are:
#       &flow_shear
#         flow_shear_moving_k = .true.
#         explicit_junction_sync = 'stage'
#         flow_shear_explicit_remap = .true.
#         moving_k_stream_defect_method = 'simple' / 'gmres' / 'auto'
#         moving_k_stream_defect_tol
#         moving_k_stream_offset_sets
#
# First the time traces of |phi|^2 and |A_parallel|^2 are compared with the
# expected output for a linear and a nonlinear simulation with the moving
# wavevector, and for the fixed-grid remap of the explicit time advance. Then
# the bookkeeping of the committed phase is checked in the restart file, a
# restarted simulation is compared with an uninterrupted one, the defect
# correction of the implicit streaming is checked to converge independently of
# its method and tolerance, a simulation with zero shear is compared with the
# fixed-grid algorithm, the remap clock is checked to stay zero without shear, and the order of the time stepping is checked by
# halving the time step.
################################################################################

# Python modules
import pytest
import shutil
import os, sys
import pathlib
import numpy as np
import xarray as xr
import netCDF4

# Package to run stella
module_path = str(pathlib.Path(__file__).parent.parent.parent / 'run_local_stella_simulation.py')
with open(module_path, 'r') as file: exec(file.read())

#-------------------------------------------------------------------------------
#                           Get the stella version                             #
#-------------------------------------------------------------------------------
@pytest.fixture(scope="session")
def stella_version(pytestconfig):
    return pytestconfig.getoption("stella_version")

#-------------------------------------------------------------------------------
#                                Helpers                                       #
#-------------------------------------------------------------------------------
def read_time_traces(netcdf_file):
    '''Return the time axis and the |phi|^2 and |A_parallel|^2 time traces.'''
    with xr.open_dataset(netcdf_file) as data:
        return np.array(data['t']), np.array(data['phi2']), np.array(data['apar2'])

def compare_with_expected_output(tmp_path, input_filename, rtol=1e-5):
    '''Compare the time traces with EXPECTED_OUTPUT.<input>.out.nc.'''
    local_netcdf_file = tmp_path / input_filename.replace('.in', '.out.nc')
    expected_netcdf_file = get_stella_expected_run_directory() / f'EXPECTED_OUTPUT.{input_filename.replace(".in","")}.out.nc'
    t1, phi1, apar1 = read_time_traces(local_netcdf_file)
    t2, phi2, apar2 = read_time_traces(expected_netcdf_file)
    assert np.allclose(t1, t2, rtol=rtol, atol=1e-20), 'The time axis does not match the expected output.'
    assert np.allclose(phi1, phi2, rtol=rtol, atol=1e-20), 'The |phi|^2 time trace does not match the expected output.'
    assert np.allclose(apar1, apar2, rtol=rtol, atol=1e-20), 'The |A_parallel|^2 time trace does not match the expected output.'
    print(f'\nThe time traces of {input_filename} match the expected output.\n')

def modified_input_file(tmp_path, input_filename, new_filename, replacements):
    '''Copy an input file to <tmp_path>/<new_filename> and replace lines in it.'''
    source = get_stella_expected_run_directory() / input_filename
    with open(source, 'r') as file: filedata = file.read()
    for old, new in replacements:
        assert filedata.count(old) == 1, f'"{old}" should occur exactly once in {input_filename}.'
        filedata = filedata.replace(old, new)
    with open(tmp_path / new_filename, 'w') as file: file.write(filedata)
    return new_filename

#-------------------------------------------------------------------------------
#               Compare the time evolution with the expected output            #
#-------------------------------------------------------------------------------
def test_whether_moving_wavevector_linear_evolves_correctly(tmp_path, stella_version):
    '''Linear electromagnetic simulation with the moving wavevector and the implicit mirror.'''
    input_filename = 'flow_shear_moving_k_linear.in'
    run_local_stella_simulation(input_filename, tmp_path, stella_version)
    compare_with_expected_output(tmp_path, input_filename)

#-------------------------------------------------------------------------------
def test_whether_moving_wavevector_nonlinear_evolves_correctly(tmp_path, stella_version):
    '''Nonlinear simulation with the moving wavevector, implicit mirror and streaming (defect correction).'''
    input_filename = 'flow_shear_moving_k_nonlinear.in'
    run_local_stella_simulation(input_filename, tmp_path, stella_version)
    compare_with_expected_output(tmp_path, input_filename)

#-------------------------------------------------------------------------------
def test_whether_explicit_remap_evolves_correctly(tmp_path, stella_version):
    '''Fixed-grid remap in a fully explicit simulation with junction synchronisation.'''
    input_filename = 'flow_shear_explicit_remap.in'
    run_local_stella_simulation(input_filename, tmp_path, stella_version)
    compare_with_expected_output(tmp_path, input_filename)

#-------------------------------------------------------------------------------
#                     Bookkeeping of the committed phase                       #
#-------------------------------------------------------------------------------
def test_whether_committed_phase_is_within_half_a_grid_spacing(tmp_path, stella_version):
    '''After nstep steps the committed phase of every row is nstep*dt wrapped into
    |tau| <= T_ky/2, with T_ky = dkx/(g_exb*ky), and the restart file carries the
    moving-k marker, the labels and the shear rate.'''
    input_filename = 'flow_shear_moving_k_restart.in'
    run_local_stella_simulation(input_filename, tmp_path, stella_version)

    # Grid, shear rate and number of steps of the input file
    with xr.open_dataset(tmp_path / input_filename.replace('.in', '.out.nc')) as data:
        akx = np.array(data['kx']); aky = np.array(data['ky'])
    g_exb = 4.0; delt = 0.002; nstep = 120

    # Read the restart file of the first processor
    restart_files = sorted((tmp_path / 'restart').glob('*.nc*'))
    assert len(restart_files) > 0, 'No restart file was written.'
    with netCDF4.Dataset(restart_files[0]) as data:
        assert int(data['moving_k'][:]) == 1, 'The moving-k marker is missing from the restart file.'
        assert int(data['moving_k_convention'][:]) == 1, 'Unexpected moving-k phase convention.'
        assert np.allclose(np.array(data['akx_label'][:]), akx, rtol=1e-12), 'The kx labels differ from the kx grid.'
        assert np.allclose(np.array(data['aky_label'][:]), aky, rtol=1e-12), 'The ky labels differ from the ky grid.'
        assert np.isclose(float(data['shear_rate'][:]), g_exb), 'The saved shear rate differs from g_exb*g_exbfac.'
        shiftstate = np.array(data['shiftstate'][:])

    # Expected committed phase: nstep*delt wrapped to [-T/2, T/2] by the relabelling
    # (the last checkpoint is written at the final step)
    dkx = np.min(np.abs(akx[np.abs(akx) > 0]))
    for iky, ky in enumerate(aky):
        if ky == 0.0:
            assert shiftstate[iky] == 0.0, 'The ky = 0 row should have zero phase.'
            continue
        period = dkx / (g_exb * ky)
        tau = nstep * delt
        tau = tau - np.floor(tau / period + 0.5) * period
        assert abs(shiftstate[iky]) <= 0.5 * period * (1 + 1e-9), f'|tau| > T/2 for ky = {ky}.'
        assert np.isclose(shiftstate[iky], tau, rtol=1e-10, atol=1e-13), f'Unexpected committed phase for ky = {ky}: {shiftstate[iky]} instead of {tau}.'
    print('\nThe committed phase and the moving-k restart variables are as expected.\n')

#-------------------------------------------------------------------------------
#                                 Restart                                      #
#-------------------------------------------------------------------------------
def test_whether_moving_wavevector_simulation_can_be_restarted(tmp_path, stella_version):
    '''A simulation of 60 + 60 steps restarted from its checkpoint (which holds a
    nonzero committed phase) must give the same time traces as 120 uninterrupted steps.'''
    input_filename = 'flow_shear_moving_k_restart.in'
    run_local_stella_simulation(input_filename, tmp_path, stella_version)
    t_full, phi_full, apar_full = read_time_traces(tmp_path / input_filename.replace('.in', '.out.nc'))

    # Run 60 steps, which writes the checkpoint, then restart up to step 120
    first = modified_input_file(tmp_path, input_filename, 'first_half.in', [('nstep = 120', 'nstep = 60')])
    run_stella(get_stella_path(stella_version), first)
    second = modified_input_file(tmp_path, input_filename, 'first_half.in',
        [("initialise_distribution_option = 'maxwellian'", "initialise_distribution_option = 'restart'")])
    run_stella(get_stella_path(stella_version), second)
    t_restarted, phi_restarted, apar_restarted = read_time_traces(tmp_path / 'first_half.out.nc')

    # Compare the last records (the restarted output file holds the second half only)
    assert np.isclose(t_restarted[-1], t_full[-1]), 'The restarted simulation did not reach the same time.'
    assert np.allclose(phi_restarted[-1], phi_full[-1], rtol=1e-8), 'The restarted |phi|^2 differs from the uninterrupted simulation.'
    assert np.allclose(apar_restarted[-1], apar_full[-1], rtol=1e-8), 'The restarted |A_parallel|^2 differs from the uninterrupted simulation.'
    print('\nThe moving-wavevector simulation was restarted correctly.\n')

#-------------------------------------------------------------------------------
#                     Defect correction of the implicit streaming              #
#-------------------------------------------------------------------------------
def test_whether_defect_correction_converges_independently_of_the_method(tmp_path, stella_version):
    '''With a tight tolerance the 'gmres' and 'auto' methods, with and without additional
    response sets, give the same solution; a loose tolerance changes it only at the level
    of the tolerance.'''
    input_filename = 'flow_shear_moving_k_linear_stream.in'
    run_local_stella_simulation(input_filename, tmp_path, stella_version)
    t_ref, phi_ref, apar_ref = read_time_traces(tmp_path / input_filename.replace('.in', '.out.nc'))
    for label, replacements, rtol in [
            ('gmres', [("moving_k_stream_defect_method = 'auto'", "moving_k_stream_defect_method = 'gmres'")], 1e-8),
            ('sets', [('moving_k_stream_offset_sets = 0', 'moving_k_stream_offset_sets = 2')], 1e-8),
            ('loose', [('moving_k_stream_defect_tol = 1e-09', 'moving_k_stream_defect_tol = 1e-5')], 1e-3)]:
        new = modified_input_file(tmp_path, input_filename, f'{label}.in', replacements)
        run_stella(get_stella_path(stella_version), new)
        t, phi, apar = read_time_traces(tmp_path / f'{label}.out.nc')
        assert np.allclose(phi, phi_ref, rtol=rtol, atol=1e-20), f'|phi|^2 with the {label} variant differs by more than {rtol}.'
        assert np.allclose(apar, apar_ref, rtol=rtol, atol=1e-20), f'|A_parallel|^2 with the {label} variant differs by more than {rtol}.'
    print('\nThe defect correction converges to the same solution for every method.\n')

#-------------------------------------------------------------------------------
#                       Zero shear reduces to the fixed-grid algorithm         #
#-------------------------------------------------------------------------------
def test_whether_zero_shear_matches_the_fixed_grid_algorithm(tmp_path, stella_version):
    '''Without shear the moving wavevector equals the fixed grid, so the moving-k
    simulation must agree with the fixed-grid remap of the explicit time advance
    (both rk4 with junction synchronisation followed by the implicit mirror split).'''
    run_local_stella_simulation('flow_shear_moving_k_zero_shear.in', tmp_path, stella_version)
    run_local_stella_simulation('flow_shear_fixed_grid_zero_shear.in', tmp_path, stella_version)
    t1, phi1, apar1 = read_time_traces(tmp_path / 'flow_shear_moving_k_zero_shear.out.nc')
    t2, phi2, apar2 = read_time_traces(tmp_path / 'flow_shear_fixed_grid_zero_shear.out.nc')
    assert np.allclose(phi1, phi2, rtol=1e-10, atol=1e-20), 'Zero-shear moving-k |phi|^2 differs from the fixed-grid algorithm.'
    assert np.allclose(apar1, apar2, rtol=1e-10, atol=1e-20), 'Zero-shear moving-k |A_parallel|^2 differs from the fixed-grid algorithm.'
    print('\nWithout shear the moving-wavevector simulation matches the fixed-grid algorithm.\n')

#-------------------------------------------------------------------------------
#                  Remap clock without shear                           #
#-------------------------------------------------------------------------------
def test_whether_remap_clock_stays_zero_without_shear(tmp_path, stella_version):
    '''With prp_shear_enabled = T and g_exb = 0 the remap is called every step but no
    mode is ever remapped, so the remap clock shift_state written to the restart file
    must stay zero; otherwise a later restart with shear would remap on its first step.'''
    input_filename = 'flow_shear_fixed_grid_zero_shear.in'
    run_local_stella_simulation(input_filename, tmp_path, stella_version)
    restart_files = sorted((tmp_path / 'restart').glob('*.nc*'))
    assert len(restart_files) > 0, 'No restart file was written.'
    with netCDF4.Dataset(restart_files[0]) as data:
        shiftstate = np.array(data['shiftstate'][:])
    assert np.all(shiftstate == 0.0), f'The remap clock advanced without shear: shiftstate = {shiftstate}.'
    print('\nThe remap clock stays zero without shear.\n')

#-------------------------------------------------------------------------------
#                          Order of the time stepping                          #
#-------------------------------------------------------------------------------
def test_whether_moving_wavevector_time_stepping_converges(tmp_path, stella_version):
    '''Halving the time step twice (same final time, fully explicit rk4 with the moving
    metric evaluated at every stage) must reduce the change of the final |phi|^2 by at
    least a factor 8. The mirror is treated explicitly here: the implicit mirror split
    is first order in time and would hide the order of the explicit scheme.'''
    input_filename = 'flow_shear_moving_k_linear.in'
    os.chdir(tmp_path)
    results = []
    for label, delt, nstep in [('dt1', 0.002, 100), ('dt2', 0.001, 200), ('dt4', 0.0005, 400)]:
        new = modified_input_file(tmp_path, input_filename, f'{label}.in',
            [('delt = 0.002', f'delt = {delt}'), ('delt_max = 0.002', f'delt_max = {delt}'), ('nstep = 200', f'nstep = {nstep}'),
             ('mirror_implicit = .true.', 'mirror_implicit = .false.')])
        run_stella(get_stella_path(stella_version), new)
        t, phi, apar = read_time_traces(tmp_path / f'{label}.out.nc')
        results.append((t[-1], phi[-1], apar[-1]))
    assert np.isclose(results[0][0], results[1][0]) and np.isclose(results[0][0], results[2][0]), 'The final times differ.'
    d1 = abs(results[0][1] - results[1][1]); d2 = abs(results[1][1] - results[2][1])
    print(f'\nChange of |phi|^2 when halving the time step: {d1:.3e} then {d2:.3e} (ratio {d1/max(d2, 1e-300):.2f}).\n')
    assert d2 < d1 / 8, 'The final |phi|^2 does not converge with the time step as expected.'
