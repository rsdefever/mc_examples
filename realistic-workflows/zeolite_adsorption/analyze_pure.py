import os
import unyt as u
import numpy as np
import matplotlib.pyplot as plt

from mosdef_cassandra.analysis import ThermoProps
from mpl_toolkits.axes_grid1.inset_locator import InsetPosition, inset_axes
from scipy.stats import linregress

plt.rc('font', family='serif')

# Make sure the range of temperatures and chemical potentials
# matches those for which the simulations were run
temperatures = [298.0 * u.K, 309.0 * u.K, 350 * u.K]
mus_adsorbate = np.arange(-49, -30, 3) * u.Unit('kJ/mol')

# Create a location to output our results
outdir="plots/"
if not os.path.isdir(outdir):
    os.mkdir(outdir)
else:
    pass

# First let's check the equilibration length. We'll plot
# the pressure and total energy for each temperature
property_ = "Energy_Total"
for temperature in temperatures:
    for mu_adsorbate in np.flip(mus_adsorbate):
        dirname = f'pure_T_{temperature:0.1f}_mu_{mu_adsorbate:.1f}'.replace(" ", "_").replace("/", "-")
        thermo_equil = ThermoProps(dirname + "/equil.out.prp")
        thermo_prod = ThermoProps(dirname + "/prod.out.prp")
        prev_plt = plt.plot(
            thermo_equil.prop("MC_STEP"),
            thermo_equil.prop(property_),
            label=mu_adsorbate
        )
        plt.plot(
            thermo_prod.prop("MC_STEP"),
            thermo_prod.prop(property_),
            color=prev_plt[0]._color
        )
    plt.title(f"T = {temperature}", fontsize=16)
    plt.xlabel("MC Step", fontsize=16, labelpad=15)
    plt.ylabel(f"{property_}, {thermo_equil.prop(property_).units}".replace("_"," "), fontsize=16, labelpad=10)
    plt.legend()
    plt.tight_layout(pad=2)

    plt.savefig(outdir + f"/{property_}_T{temperature.value}K.pdf")
    plt.close()

property_ = "Pressure"
for temperature in temperatures:
    for mu_adsorbate in np.flip(mus_adsorbate):
        dirname = f'pure_T_{temperature:0.1f}_mu_{mu_adsorbate:.1f}'.replace(" ", "_").replace("/", "-")
        thermo_equil = ThermoProps(dirname + "/equil.out.prp")
        thermo_prod = ThermoProps(dirname + "/prod.out.prp")
        prev_plt = plt.plot(
            thermo_equil.prop("MC_STEP"),
            thermo_equil.prop(property_),
            label=mu_adsorbate
        )
        plt.plot(
            thermo_prod.prop("MC_STEP"),
            thermo_prod.prop(property_),
            color=prev_plt[0]._color
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
    for mu_adsorbate in np.flip(mus_adsorbate):
        dirname = f'pure_T_{temperature:0.1f}_mu_{mu_adsorbate:.1f}'.replace(" ", "_").replace("/", "-")
        thermo = ThermoProps(dirname + "/prod.out.prp")
        width=5
        n = np.ceil((thermo.prop("Nmols").value.max() - thermo.prop("Nmols").value.min())/width)
        parent_ax.hist(
            thermo.prop("Nmols"),
            alpha=0.6,
            align="left",
            bins=int(n),
            label=mu_adsorbate,
            density=True
        )
        width=1
        n = np.ceil((thermo.prop("Nmols").value.max() - thermo.prop("Nmols").value.min())/width)
        inset_ax.hist(
            thermo.prop("Nmols"),
            alpha=0.6,
            align="left",
            bins=int(n),
            label=mu_adsorbate,
            density=True
        )
   
    #parent_ax.set_title(f"T = {temperature:.0f}", fontsize=16)
    parent_ax.set_xlabel("Number of molecules", fontsize=16, labelpad=10)
    parent_ax.set_ylabel("Probability", fontsize=16, labelpad=10)
    parent_ax.tick_params(axis='both', labelsize=12)

    inset_ax.set_xlim(-2,12)
    inset_ax.set_xlabel("Number of molecules")
    inset_ax.set_ylabel("Probability")
    parent_ax.legend(loc="upper left")
    fig.tight_layout(pad=2)
    fig.savefig(outdir + f"/nmols-hist_T{temperature.value}K.pdf")
    plt.close()

# Finally let's plot the relationship between the chemical potential
# and pressure at both temperatures
for temperature in temperatures:
    pressures = []
    for mu_adsorbate in mus_adsorbate:
        dirname = f'pure_T_{temperature:0.1f}_mu_{mu_adsorbate:.1f}'.replace(" ", "_").replace("/", "-")
        thermo = ThermoProps(dirname + "/prod.out.prp")
        pressures.append(np.mean(thermo.prop("Pressure", start=100000)))
    pressures = u.unyt_array(pressures)

    # Fit a line to mu vs. P
    slope, intercept, r_value, p_value, stderr = linregress(
        mus_adsorbate.to_value('kJ/mol').flatten(),
        y=np.log(pressures.to_value(u.bar)).flatten()
    )
    # Determine chemical potentials
    fit_press = np.exp(slope * np.arange(-50,-29,1) + intercept) * u.bar

    prev_plt = plt.plot(
        mus_adsorbate.to_value('kJ/mol'),
        pressures.to_value('bar'), 'o',
        markersize=10,
    )
    plt.plot(
        np.arange(-50,-29,1),
        fit_press.to_value('bar'),
        '-',
        color=prev_plt[0]._color,
        label=f"{temperature:.0f}"
    )
    plt.xlabel("Chemical potential (kJ/mol)", fontsize=16, labelpad=15)
    plt.ylabel("Pressure (bar)", fontsize=16, labelpad=10)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.yscale('log')
    plt.xlim(-50,-30)
    plt.legend(fontsize=16)
    plt.tight_layout(pad=2.0)

    plt.savefig(outdir + f"/chempot_vs_pressure.pdf")

