################################################################################
#          Check that an embedded coarse kx grid stays on its comb             #
################################################################################
# A linear state of a jtwist = 2 grid is embedded in a jtwist = 4 grid by
# keeping only the even kx labels (the odd labels are set to zero in a
# checkpoint). In a linear simulation every operation acts within a ky row,
# shifts a whole row by one label (relabelling with the moving wavevector, or
# the fixed-grid remap) or couples labels four grid spacings apart (the
# twist-and-shift chains of this grid), so the energy on the labels off this
# comb must stay at roundoff for both treatments of the flow shear. The run
# covers about two relabel periods of the lowest ky.
################################################################################

# Python modules
import pytest
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
def modified_input_file(tmp_path, input_filename, new_filename, replacements):
    '''Copy an input file to <tmp_path>/<new_filename> and replace lines in it.'''
    source = get_stella_expected_run_directory() / input_filename
    with open(source, 'r') as file: filedata = file.read()
    for old, new in replacements:
        assert filedata.count(old) == 1, f'"{old}" should occur exactly once in {input_filename}.'
        filedata = filedata.replace(old, new)
    with open(tmp_path / new_filename, 'w') as file: file.write(filedata)
    return new_filename

def keep_even_kx_labels(restart_file, akx, naky):
    '''Set g to zero for the odd kx labels in a single-processor checkpoint
    (layout yxzs: the kx index of global point i is 1 + (i // naky) mod nakx).'''
    nakx = len(akx)
    dkx = np.min(np.abs(akx[np.abs(akx) > 0]))
    even = (np.rint(akx / dkx).astype(int) % 2 == 0)
    with netCDF4.Dataset(restart_file, 'r+') as data:
        nglo = data['gr'].shape[0]
        ikx = (np.arange(nglo) // naky) % nakx
        odd = ~even[ikx]
        for name in ('gr', 'gi'):
            g = data[name][:]
            g[odd, :, :] = 0.0
            data[name][:] = g
    return even

def comb_energy_ratio(netcdf_file, even):
    '''Largest, over time and non-zonal ky, of the |phi|^2 energy of the parity
    class with less energy divided by that of the class with more energy. Every
    relabelling or remap flips the parity of the comb of a row, so the class with
    less energy is the one off the comb. Rows whose labels have all left the kx
    grid (large ky, nothing enters in a linear simulation) are skipped once empty.'''
    with xr.open_dataset(netcdf_file) as data:
        phi = np.array(data['phi_vs_t'])          # (t, tube, z, kx, ky, ri)
        ky = np.array(data['ky'])
    energy = np.sum(phi[..., 0]**2 + phi[..., 1]**2, axis=(1, 2))   # (t, kx, ky)
    worst = 0.0
    for iky in range(len(ky)):
        if ky[iky] == 0.0: continue
        e_even = np.sum(energy[:, even, iky], axis=1)
        e_odd = np.sum(energy[:, ~even, iky], axis=1)
        on = np.maximum(e_even, e_odd); off = np.minimum(e_even, e_odd)
        if iky == 1: assert np.all(on > 0.0), f'No energy left in the ky = {ky[iky]} row.'
        filled = on > 0.0
        worst = max(worst, float(np.max(off[filled] / on[filled])))
    return worst

#-------------------------------------------------------------------------------
#                       Comb invariance of the linear dynamics                 #
#-------------------------------------------------------------------------------
def test_whether_embedded_coarse_grid_stays_on_its_comb(tmp_path, stella_version):
    '''Moving wavevector and fixed-grid remap, restarted from the same checkpoint
    holding only the even kx labels; the off-comb energy must stay at roundoff.'''
    input_filename = 'flow_shear_comb_jtwist4.in'
    os.chdir(tmp_path)

    # Checkpoint after one step without shear (zero remap phase), one processor
    seed = modified_input_file(tmp_path, input_filename, 'seed.in', [('nstep = 160', 'nstep = 1')])
    run_stella(get_stella_path(stella_version), seed, nproc=1)
    with xr.open_dataset(tmp_path / 'seed.out.nc') as data:
        akx = np.array(data['kx']); naky = len(data['ky'])
    restart_files = sorted((tmp_path / 'restart').glob('*.nc*'))
    assert len(restart_files) == 1, 'Expected a single restart file.'
    even = keep_even_kx_labels(restart_files[0], akx, naky)

    # Restart with shear: moving wavevector, and fixed-grid remap of the explicit time advance
    common = [("initialise_distribution_option = 'maxwellian'", "initialise_distribution_option = 'restart'"),
              ("restart_dir = 'restart'", "restart_dir = 'restart'\n    restart_file = 'seed.nc'"),
              ('g_exb = 0.0', 'g_exb = 4.0'), ('save_for_restart = .true.', 'save_for_restart = .false.')]
    moving = modified_input_file(tmp_path, input_filename, 'moving_k.in', common
        + [('flow_shear_moving_k = .false.', 'flow_shear_moving_k = .true.'),
           ('flow_shear_explicit_remap = .true.', 'flow_shear_moving_k_import_zero_phase_seed = .true.')])
    remap = modified_input_file(tmp_path, input_filename, 'remap.in', common)
    for name in (moving, remap):
        run_stella(get_stella_path(stella_version), name, nproc=1)
        ratio = comb_energy_ratio(tmp_path / name.replace('.in', '.out.nc'), even)
        print(f'\n{name}: largest off-comb to on-comb energy ratio = {ratio:.3e}\n')
        assert ratio < 1e-12, f'{name}: the energy leaves the comb of the embedded coarse grid (ratio {ratio:.3e}).'
