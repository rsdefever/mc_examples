import os
import unyt as u
import numpy as np
import matplotlib.pyplot as plt

from mosdef_cassandra.analysis import ThermoProps
from mpl_toolkits.axes_grid1.inset_locator import InsetPosition, inset_axes

# Make sure the range of temperatures and chemical potentials
# matches those for which the simulations were run
temperatures = [298.0 * u.K, 309.0 * u.K]
mus_adsorbate = np.arange(-46, -25, 3) * u.Unit('kJ/mol')

# Create a location to output our results
outdir="plots/"
if not os.path.isdir(outdir):
    os.mkdir(outdir)
else:
    pass

# First let's check the equilibration length. We'll plot
# the pressure and total energy for each temperature
for property_ in ["Energy_Total", "Pressure"]:
    for temperature in temperatures:
        for mu_adsorbate in mus_adsorbate:
            dirname = f'pure_T_{temperature:0.1f}_mu_{mu_adsorbate:.1f}'.replace(" ", "_").replace("/", "-")
            thermo_equil = ThermoProps(dirname + "/equil.out.prp")
            thermo_prod = ThermoProps(dirname + "/prod.out.prp")
            plt.plot(
                thermo_equil.prop("MC_STEP"),
                thermo_equil.prop(property_),
                label=mu_adsorbate
            )
            plt.plot(
                thermo_prod.prop("MC_STEP"),
                thermo_prod.prop(property_),
            )
        plt.title(f"{property_}, T = {temperature}")
        plt.xlabel("MC Step")
        plt.ylabel(f"{property_}, {thermo_equil.prop(property_).units}")
        plt.legend()

        plt.savefig(outdir + f"/{property_}_T{temperature.value}K.pdf")
        plt.close()

# Next let's look at the number of molecules in the
# simulation boxes to make sure our boxes were large enough.
# We will only look at the production simulations
for temperature in temperatures:
    fig, parent_ax = plt.subplots()
    inset_ax = inset_axes(
        parent_ax, width="45%", height=1.0, loc="upper right",
        bbox_to_anchor=(-0.02,-0.02,1,1), bbox_transform=parent_ax.transAxes
    )
    for mu_adsorbate in mus_adsorbate:
        dirname = f'pure_T_{temperature:0.1f}_mu_{mu_adsorbate:.1f}'.replace(" ", "_").replace("/", "-")
        thermo = ThermoProps(dirname + "/prod.out.prp")
        parent_ax.hist(thermo.prop("Nmols"), alpha=0.75, align="left")
        inset_ax.hist(thermo.prop("Nmols"), alpha=0.75, align="left")
    
    parent_ax.set_title("Nmols in simulation box")
    parent_ax.set_xlabel("Number of molecules")
    parent_ax.set_ylabel("Frequency")

    inset_ax.set_xlim(-2,10)
    inset_ax.set_xlabel("Number of molecules")
    inset_ax.set_ylabel("Frequency")
    fig.savefig(outdir + f"/nmols-hist_T{temperature.value}K.pdf")

# Finally let's plot the relationship between the chemical potential
# and pressure at both temperatures
for temperature in temperatures:
    pressures = []
    for mu_adsorbate in mus_adsorbate:
        dirname = f'pure_T_{temperature:0.1f}_mu_{mu_adsorbate:.1f}'.replace(" ", "_").replace("/", "-")
        thermo = ThermoProps(dirname + "/prod.out.prp")
        pressures.append(np.mean(thermo.prop("Pressure")))
    pressures = u.unyt_array(pressures)
    plt.plot(
        mus_adsorbate.to_value('kJ/mol'),
        pressures.to_value('bar'), 'o-'
    )
    plt.xlabel("Chemical potential (kJ/mol)")
    plt.ylabel("Pressure (bar)")
    plt.yscale('log')

    plt.savefig(outdir + f"/chempot_vs_pressure.pdf")

