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


def main():

    # Load TON from a CIF file, replicate the cell
    # Use mbuild to create a zeolite supercell from CIF
    lattice = mbuild.lattice.load_cif("zeolite_resources/structures/TON.cif")
    compound_dict = {
        "Si": mbuild.Compound(name="Si"),
        "O": mbuild.Compound(name="O"),
    }
    zeolite = lattice.populate(compound_dict, 2, 2, 6)
    
    # Create a CG methane, load and apply ff
    methane = mbuild.Compound(name="_CH4")
    ff_ads = foyer.Forcefield("zeolite_resources/ffxml/adsorbates.xml")
    methane_ff = ff_ads.apply(methane)

    # Define pure fluid temperatures and chemical potentials
    temperatures = [298 * u.K, 309 * u.K, 350 * u.K]
    mus_fluid = np.arange(-49, -30, 3) * u.Unit("kJ/mol")

    # Define the pressures at which we wish to study adsorption
    pressures = [
        0.01,
        0.1,
        0.25,
        0.5,
        0.75,
        1.0,
        2.0,
        3.0,
        5.0,
    ] * u.bar

    # Select the zeolite ff
    zeo_ff_names = ["june", "trappe"]

    # Define a few custom_args that will be
    # the same for all zeolite simulations
    custom_args = {
        "charge_style": "none",
        "vdw_cutoff": 14.0 * u.angstrom,
        "prop_freq": 10,
        "max_molecules": [1, 10000],
    }

    # Loop over different zeolite ff's
    for zeo_ff_name in zeo_ff_names:

        # Load and apply ff to the zeolite structure
        ff_zeo = foyer.Forcefield(f"zeolite_resources/ffxml/zeo_{zeo_ff_name}.xml")
        zeolite_ff = ff_zeo.apply(zeolite)

        # Create the box_list, species_list, System, and MoveSet.
        # These are not dependent upon (T,P) condition
        box_list = [zeolite]
        species_list = [zeolite_ff, methane_ff]
        mols_in_boxes = [[1, 0]]

        system = mc.System(box_list, species_list, mols_in_boxes=mols_in_boxes)
        moveset = mc.MoveSet("gcmc", species_list)

        # Loop over each temperature to compute an isotherm
        for temperature in temperatures:

            # Before we begin we must determine the
            # chemical potentials required to achieve
            # the desired pressures
            fluid_pressures = []
            for mu_fluid in mus_fluid:
                dirname = f"fluid_T_{temperature:0.1f}_mu_{mu_fluid:.1f}".replace(
                    " ", "_"
                ).replace("/", "-")
                thermo = ThermoProps(dirname + "/prod.out.prp")
                fluid_pressures.append(np.mean(thermo.prop("Pressure")))
            fluid_pressures = u.unyt_array(fluid_pressures)

            # Fit a line to mu vs. P
            slope, intercept, r_value, p_value, stderr = linregress(
                np.log(fluid_pressures.to_value(u.bar)).flatten(),
                y=mus_fluid.to_value("kJ/mol").flatten(),
            )
            # Determine chemical potentials
            mus = (slope * np.log(pressures.in_units(u.bar)) + intercept) * u.Unit(
                "kJ/mol"
            )

            # Loop over each pressure and run the MC simulation!
            for (pressure, mu) in zip(pressures, mus):
                print(f"\nRun simulation: T = {temperature}, P = {pressure}\n")
                dirname = f"zeo_ff_{zeo_ff_name}_T_{temperature:0.1f}_P_{pressure:0.2f}".replace(
                    " ", "_"
                ).replace(
                    "/", "-"
                )
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
                        chemical_potentials=["none", mu],
                        **custom_args,
                    )

                    mc.restart(
                        system=system,
                        moveset=moveset,
                        run_type="prod",
                        run_length=200000,
                        temperature=temperature,
                        run_name="prod",
                        restart_name="equil",
                        chemical_potentials=["none", mu],
                        **custom_args,
                    )


if __name__ == "__main__":
    main()

