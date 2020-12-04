import os
import mbuild
import foyer
import mosdef_cassandra as mc
import unyt as u
import numpy as np
import matplotlib.pyplot as plt

from scipy.stats import linregress
from mosdef_cassandra.analysis import ThermoProps
from mosdef_cassandra.utils.tempdir import temporary_cd
from pkg_resources import resource_filename


def main():
    # Create a CG methane, load and apply ff
    methane = mbuild.Compound(name="_CH4")
    ff_path = resource_filename(
        "mc_examples",
        "realistic_workflows/zeolite_adsorption/resources/ffxml/adsorbates.xml",
    )
    ff_ads = foyer.Forcefield(ff_path)
    methane_ff = ff_ads.apply(methane)

    # Define a few keyword args that will be the
    # same for all fluid phase simulations
    custom_args = {
        "charge_style": "none",
        "vdw_cutoff": 14.0 * u.angstrom,
        "prop_freq": 10,
    }

    # Define temperatures and chemical potential values
    # for the single phase GCMC simulations
    temperatures = [298 * u.K, 309 * u.K, 350 * u.K]
    mus = np.arange(-49, -30, 3) * u.Unit("kJ/mol")

    # Loop over temperatures and mus and run simulations
    for temperature in temperatures:
        for mu in mus:

            # For each simulation -- first estimate
            # the box volume to target ~40 molecules
            n_methane_target = 40
            beta = 1.0 / (u.kb * temperature)
            mass = 16.04 * u.amu
            debroglie = np.sqrt(2 * np.pi * u.hbar ** 2 * beta / mass)
            vol = n_methane_target * debroglie ** 3 * np.exp(-beta.to("mol/kJ") * mu)
            boxl = (vol ** (1.0 / 3.0)).to_value("nm")
            if boxl < 2.0 * custom_args["vdw_cutoff"].to_value("nm"):
                boxl = 2.0 * custom_args["vdw_cutoff"].to_value("nm")

            # Define the species list, box list, system, moveset
            # We start with an empty box
            species_list = [methane_ff]
            box_list = [mbuild.Box([boxl, boxl, boxl])]
            system = mc.System(box_list, species_list)
            moveset = mc.MoveSet("gcmc", species_list)

            # Create a new directory, temporary cd, and run
            # the equilibration and production simulations!
            print(f"\nRun simulation: T = {temperature}, mu = {mu}\n")
            dirname = f"fluid_T_{temperature:0.1f}_mu_{mu:.1f}".replace(
                " ", "_"
            ).replace("/", "-")
            if not os.path.isdir(dirname):
                os.mkdir(dirname)
            else:
                pass
            with temporary_cd(dirname):
                mc.run(
                    system=system,
                    moveset=moveset,
                    run_type="equil",
                    run_length=50000,
                    temperature=temperature,
                    run_name="equil",
                    chemical_potentials=[mu],
                    **custom_args,
                )

                mc.restart(
                    system=system,
                    moveset=moveset,
                    run_type="prod",
                    run_length=400000,
                    temperature=temperature,
                    run_name="prod",
                    restart_name="equil",
                    chemical_potentials=[mu],
                    **custom_args,
                )


if __name__ == "__main__":
    main()
