import os
import unyt as u
import numpy as np
import matplotlib.pyplot as plt
import seaborn

from mosdef_cassandra.analysis import ThermoProps
from scipy.stats import linregress
from matplotlib.ticker import MultipleLocator, AutoMinorLocator
from matplotlib import rcParams

rcParams["font.sans-serif"] = "Arial"
rcParams["font.family"] = "sans-serif"


def main():

    # Make sure the range of temperatures and chemical potentials
    # matches those for which the simulations were run
    temperatures = [298 * u.K, 309 * u.K, 350 * u.K]
    mus = np.arange(-49, -30, 3) * u.Unit("kJ/mol")

    # Create a location to output our results
    outdir = "plots/"
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    else:
        pass

    # First let's check the equilibration length. We'll plot
    # the pressure and total energy for each temperature
    for property_ in ["Energy_Total", "Pressure"]:
        for temperature in temperatures:
            fig, ax = plt.subplots()
            for mu in np.flip(mus):
                dirname = f"fluid_T_{temperature:0.1f}_mu_{mu:.1f}".replace(
                    " ", "_"
                ).replace("/", "-")
                thermo_equil = ThermoProps(dirname + "/equil.out.prp")
                thermo_prod = ThermoProps(dirname + "/prod.out.prp")
                prev_plt = ax.plot(
                    thermo_equil.prop("MC_STEP"),
                    thermo_equil.prop(property_),
                    label=mu,
                )
                ax.plot(
                    thermo_prod.prop("MC_STEP"),
                    thermo_prod.prop(property_),
                    color=prev_plt[0]._color,
                )
            ax.set_title(f"T = {temperature}", fontsize=16)
            ax.set_xlabel("MC Step", fontsize=16, labelpad=15)
            ax.set_ylabel(
                f"{property_}, {thermo_equil.prop(property_).units}".replace("_", " "),
                fontsize=16,
                labelpad=10,
            )
            ax.tick_params(
                which="both", direction="in", labelsize=14, top=True, right=True
            )
            fig.legend()
            fig.tight_layout(pad=2)
            fig.savefig(outdir + f"/{property_}_T{temperature.value}K.pdf")

    # Next let's look at the number of molecules in the
    # simulation boxes to make sure our boxes were large enough.
    # We will only look at the production simulations
    for temperature in temperatures:
        fig, ax = plt.subplots()
        for mu in np.flip(mus):
            dirname = f"fluid_T_{temperature:0.1f}_mu_{mu:.1f}".replace(
                " ", "_"
            ).replace("/", "-")
            thermo = ThermoProps(dirname + "/prod.out.prp")
            width = 5
            ax.hist(
                thermo.prop("Nmols"),
                alpha=0.4,
                align="left",
                bins=range(15, 70, 3),
                label=str(mu),
                density=True,
            )

        ax.set_xlabel("Number of molecules", fontsize=16, labelpad=10)
        ax.set_ylabel("Probability", fontsize=16, labelpad=10)
        ax.tick_params(which="both", direction="in", labelsize=14, top=True, right=True)
        ax.xaxis.set_minor_locator(AutoMinorLocator(n=2))
        ax.yaxis.set_minor_locator(AutoMinorLocator(n=2))
        fig.tight_layout(pad=2)
        fig.legend(fontsize=16, loc="upper right", bbox_to_anchor=(0.95, 0.95))
        fig.savefig(outdir + f"/nmols-hist_T{temperature.value}K.pdf")

    # Finally let's plot the relationship between mu
    # and pressure at all temperatures
    fig, ax = plt.subplots()
    for temperature in temperatures:
        pressures = []
        for mu in mus:
            dirname = f"fluid_T_{temperature:0.1f}_mu_{mu:.1f}".replace(
                " ", "_"
            ).replace("/", "-")
            thermo = ThermoProps(dirname + "/prod.out.prp")
            pressures.append(np.mean(thermo.prop("Pressure", start=100000)))
        pressures = u.unyt_array(pressures)

        # Fit a line to mu vs. P
        slope, intercept, r_value, p_value, stderr = linregress(
            mus.to_value("kJ/mol").flatten(),
            y=np.log(pressures.to_value(u.bar)).flatten(),
        )
        # Determine chemical potentials
        fit_press = np.exp(slope * np.arange(-50, -29, 1) + intercept) * u.bar

        # Plot
        prev_plt = ax.plot(
            mus.to_value("kJ/mol"),
            pressures.to_value("bar"),
            "o",
            markersize=10,
        )
        ax.plot(
            np.arange(-50, -29, 1),
            fit_press.to_value("bar"),
            "-",
            color=prev_plt[0]._color,
            label=f"{temperature:.0f}",
        )
    ax.set_xlabel("Chemical potential (kJ/mol)", fontsize=16, labelpad=15)
    ax.xaxis.set_minor_locator(AutoMinorLocator(3))
    ax.xaxis.set_major_locator(MultipleLocator(3))
    ax.set_ylabel("Pressure (bar)", fontsize=16, labelpad=10)
    ax.tick_params(which="both", direction="in", labelsize=14, top=True, right=True)
    ax.set_yscale("log")
    ax.set_xlim(-50, -30)
    fig.legend(fontsize=16, loc="upper left", bbox_to_anchor=(0.15, 0.95))
    fig.tight_layout(pad=2.0)
    fig.savefig(outdir + f"/chempot_vs_pressure.pdf")


if __name__ == "__main__":
    main()
