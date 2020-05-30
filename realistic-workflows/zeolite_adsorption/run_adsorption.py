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

# Define the same temperatures and the pure phase chemical potentials
temperatures = [298.0 * u.K, 309.0 * u.K, 350 * u.K]
mus_adsorbate = np.arange(-49, -30, 3) * u.Unit('kJ/mol')

# Define the pressures at which we wish to study adsorption
pressures = [
    0.01,
    0.1,
    0.25,
    0.5,
    0.75
    1.0,
    2.0,
    3.0,
    5.0,
] * u.bar

# Create a CG methane, load and apply ff
methane = mbuild.Compound(name='_CH4')
ff_ads = foyer.Forcefield('resources/ffxml/adsorbates.xml')
methane_ff = ff_ads.apply(methane)

# Load the zeolite, load and apply ff
zeolite = mbuild.load("resources/structures/TON_2x2x6.pdb")
ff_zeo = foyer.Forcefield("resources/ffxml/zeo_june.xml")
zeolite_ff = ff_zeo.apply(zeolite)

# Define a few custom_args that will be the same for
# all zeolite simulations
custom_args = {
    "charge_style" : "none",
    "vdw_cutoff" : 14.0 * u.angstrom,
    "prop_freq" : 10,
    "max_molecules" : [1, 10000]
}

# Create the box_list, species_list, System, and MoveSet
# since they will be the same for each condition
box_list = [zeolite]
species_list = [zeolite_ff, methane_ff]
mols_in_boxes = [[1,0]]

system = mc.System(box_list, species_list, mols_in_boxes=mols_in_boxes)
moveset = mc.MoveSet('gcmc', species_list)

# Let the fun begin
for temperature in temperatures:
    # Before we begin we must do a bit more analysis to determine the
    # chemical potentials required to achieve the desired pressures
    pure_pressures = []
    for mu_adsorbate in mus_adsorbate:
        dirname = f'pure_T_{temperature:0.1f}_mu_{mu_adsorbate:.1f}'.replace(" ", "_").replace("/", "-")
        thermo = ThermoProps(dirname + "/prod.out.prp")
        pure_pressures.append(np.mean(thermo.prop("Pressure")))
    pure_pressures = u.unyt_array(pure_pressures)

    # Fit a line to mu vs. P
    slope, intercept, r_value, p_value, stderr = linregress(
        np.log(pure_pressures.to_value(u.bar)).flatten(),
        y=mus_adsorbate.to_value('kJ/mol').flatten()
    )
    # Determine chemical potentials
    mus = (slope * np.log(pressures.in_units(u.bar)) + intercept) * u.Unit('kJ/mol')

    # Now loop over each pressure and run the MC simulation!
    for (pressure, mu) in zip(pressures, mus):
        print(f"\nRun simulation: T = {temperature}, P = {pressure}\n")
        dirname = f'zeo_T_{temperature:0.1f}_P_{pressure:0.2f}'.replace(" ", "_").replace("/", "-")
        if not os.path.isdir(dirname):
            os.mkdir(dirname)
        else:
            pass
        with temporary_cd(dirname):

            mc.run(
                system=system,
                moveset=moveset,
                run_type="equil",
                run_length=30000,
                temperature=temperature,
                run_name='equil',
                chemical_potentials = ["none", mu],
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
                chemical_potentials = ["none", mu],
                **custom_args,
            )
