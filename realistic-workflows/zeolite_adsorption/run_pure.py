import mbuild
import foyer
import mosdef_cassandra as mc
import unyt as u

import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

from mosdef_cassandra.analysis import ThermoProps
from mosdef_cassandra.utils.tempdir import temporary_cd

# Filter some warnings -- to cleanup output for this demo
from warnings import filterwarnings
filterwarnings('ignore', category=UserWarning)
from parmed.exceptions import OpenMMWarning
filterwarnings('ignore', category=OpenMMWarning)

# Create a CG methane, load and apply ff
methane = mbuild.Compound(name='_CH4')
ff_ads = foyer.Forcefield('resources/ffxml/adsorbates.xml')
methane_typed = ff_ads.apply(methane)

# Define a few custom_args that will be the same for
# all pure phase simulations
custom_args = {
    "charge_style" : "none",
    "vdw_cutoff" : 14.0 * u.angstrom,
    "prop_freq" : 10,
    "max_molecules" : [20000]
}

# Run pure phase simulations at two temperatures
# Use the same range of chemical potentials at both temperatures
temperatures = [298.0 * u.K, 309.0 * u.K, 350 * u.K]
mus_adsorbate = np.arange(-49, -30, 3) * u.Unit('kJ/mol')

for temperature in temperatures:
    for mu_adsorbate in mus_adsorbate:
        print(f"\nRun simulation: T = {temperature}, mu = {mu_adsorbate}\n")
        dirname = f'pure_T_{temperature:0.1f}_mu_{mu_adsorbate:.1f}'.replace(" ", "_").replace("/", "-")
        if not os.path.isdir(dirname):
            os.mkdir(dirname)
        else:
            pass
        with temporary_cd(dirname):
            species_list = [methane_typed]
            if mu_adsorbate < -34:
                boxl = 22. # nm
            else:
                boxl = 5. # nm
            box_list = [mbuild.Box([boxl,boxl,boxl])]
            system = mc.System(box_list, species_list)
            moveset = mc.MoveSet('gcmc', species_list)
    
            mc.run(
                system=system,
                moveset=moveset,
                run_type="equil",
                run_length=30000,
                temperature=temperature,
                run_name='equil',
                chemical_potentials = [mu_adsorbate],
                **custom_args,
            )
    
            mc.restart(
                system=system,
                moveset=moveset,
                run_type="prod",
                run_length=100000,
                temperature=temperature,
                run_name='prod',
                restart_name='equil',
                chemical_potentials = [mu_adsorbate],
                **custom_args,
            )
