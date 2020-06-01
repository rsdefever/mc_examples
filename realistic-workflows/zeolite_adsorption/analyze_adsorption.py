import os
import unyt as u
import numpy as np
import matplotlib.pyplot as plt

from mosdef_cassandra.analysis import ThermoProps
from mpl_toolkits.axes_grid1.inset_locator import InsetPosition, inset_axes

# Make sure the range of temperatures and pressures
# matches those for which the simulations were run
temperatures = [298.0 * u.K, 309.0 * u.K, 350 * u.K]
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

# Define the number of unit cells (2x2x6=24)
n_unitcells = 24

# Create a location to output our results
outdir="plots/"
if not os.path.isdir(outdir):
    os.mkdir(outdir)
else:
    pass

# First let's check the equilibration length. We'll plot
# the number of molecules/unitcell over the simulations

property_ = "Nmols_2"
for temperature in temperatures:
    for pressure in pressures:
        dirname = f'zeo_T_{temperature:0.1f}_P_{pressure:0.2f}'.replace(" ", "_").replace("/", "-")
        thermo_equil = ThermoProps(dirname + "/equil.out.prp")
        thermo_prod = ThermoProps(dirname + "/prod.out.prp")
        plt.plot(
            np.hstack((thermo_equil.prop("MC_STEP"),thermo_prod.prop("MC_STEP"))),
            np.hstack((thermo_equil.prop(property_),thermo_prod.prop(property_)))/n_unitcells,
            label=pressure
        )
    plt.title("$\mathregular{N_{methane}}$" + f", T = {temperature}")
    plt.xlabel("MC Step")
    plt.ylabel("$\mathregular{N_{methane}}$/unit cell")
    plt.legend()

    plt.savefig(outdir + f"/Nmols_T{temperature.value}K.pdf")
    plt.close()

# Next let's calculate the average number of molecules per unit
# cell during the production simulation and plot the isotherms

# Load lit results to compare
lit_298K = np.genfromtxt("resources/lit_results/tjune_TON-methane_298K.txt", skip_header=1)
lit_309K = np.genfromtxt("resources/lit_results/tjune_TON-methane_309K.txt", skip_header=1)

for temperature in temperatures:
    avg_n = []
    for pressure in pressures:
        dirname = f'zeo_T_{temperature:0.1f}_P_{pressure:0.2f}'.replace(" ", "_").replace("/", "-")
        thermo = ThermoProps(dirname + "/prod.out.prp")
        avg_n.append(np.mean(thermo.prop("Nmols_2"))/n_unitcells)
    avg_n = u.unyt_array(avg_n)
    plt.scatter(
        pressures.to_value('bar'),
        avg_n.value,
        label=temperature,
        s=100,
    )

plt.scatter(
    lit_298K[:,0],
    lit_298K[:,1],
    marker='1',
    label="Romanielo 2020 298.0 K",
    s=150
)
plt.scatter(
    lit_309K[:,0],
    lit_309K[:,1],
    marker='2',
    label="Romanielo 2020 309.0 K",
    s=150
)

plt.ylabel("$\mathregular{N_{methane}}$/unit cell", fontsize=24, labelpad=10)
plt.xlabel("Pressure (bar)", fontsize=24, labelpad=10)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.xscale('log')
plt.legend(fontsize=16)
plt.tight_layout()

plt.savefig(outdir + f"/isotherm.pdf")


