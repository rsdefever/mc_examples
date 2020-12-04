import os
import unyt as u
import numpy as np
import matplotlib.pyplot as plt
import seaborn

from mosdef_cassandra.analysis import ThermoProps
from scipy.stats import linregress
from matplotlib.ticker import MultipleLocator, AutoMinorLocator
from matplotlib import rcParams
from pkg_resources import resource_filename

rcParams["font.sans-serif"] = "Arial"
rcParams["font.family"] = "sans-serif"


def main():

    # Make sure the range of temperatures and pressures
    # matches those for which the simulations were run
    zeo_ff_names = ["june", "trappe"]
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
    outdir = "plots/"
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    else:
        pass

    # First let's check the equilibration length. We'll plot
    # the number of molecules/unitcell over the simulations
    for zeo_ff_name in zeo_ff_names:
        property_ = "Nmols_2"
        for temperature in temperatures:
            fig, ax = plt.subplots()
            for pressure in pressures:
                dirname = f"zeo_ff_{zeo_ff_name}_T_{temperature:0.1f}_P_{pressure:0.2f}".replace(
                    " ", "_"
                ).replace(
                    "/", "-"
                )
                thermo_equil = ThermoProps(dirname + "/equil.out.prp")
                thermo_prod = ThermoProps(dirname + "/prod.out.prp")
                ax.plot(
                    np.hstack(
                        (thermo_equil.prop("MC_STEP"), thermo_prod.prop("MC_STEP"))
                    ),
                    np.hstack(
                        (thermo_equil.prop(property_), thermo_prod.prop(property_))
                    )
                    / n_unitcells,
                    label=pressure,
                )
            ax.set_title(
                "$\mathregular{N_{methane}}$" + f", T = {temperature}", fontsize=16
            )
            ax.set_xlabel("MC Step", fontsize=16)
            ax.set_ylabel("$\mathregular{N_{methane}}$/unit cell", fontsize=16)
            ax.tick_params(
                which="both", direction="in", labelsize=14, top=True, right=True
            )
            fig.legend()
            fig.tight_layout(pad=2)
            fig.savefig(
                outdir + f"/Nmols_zeo_ff_{zeo_ff_name}_T{temperature.value}K.pdf"
            )

    # Next let's calculate the average number of molecules per unit
    # cell during the production simulation and plot the isotherms
    fig, ax = plt.subplots()

    # Load lit results to compare
    file_path = resource_filename(
        "mc_examples",
        "realistic_workflows/zeolite_adsorption/resources/lit_results/tjune_TON-methane_298K.txt",
    )
    lit_298K = np.genfromtxt(file_path, skip_header=1)
    file_path = resource_filename(
        "mc_examples",
        "realistic_workflows/zeolite_adsorption/resources/lit_results/tjune_TON-methane_309K.txt",
    )
    lit_309K = np.genfromtxt(file_path, skip_header=1)

    for temperature in temperatures:
        zeo_ff_name = "june"
        avg_n = []
        for pressure in pressures:
            dirname = (
                f"zeo_ff_{zeo_ff_name}_T_{temperature:0.1f}_P_{pressure:0.2f}".replace(
                    " ", "_"
                ).replace("/", "-")
            )
            thermo = ThermoProps(dirname + "/prod.out.prp")
            avg_n.append(np.mean(thermo.prop("Nmols_2")) / n_unitcells)
        avg_n = u.unyt_array(avg_n)
        prev_plt = ax.scatter(
            pressures.to_value("bar"),
            avg_n.value,
            label=f"{temperature:.0f} June",
            s=100,
        )

        zeo_ff_name = "trappe"
        avg_n = []
        for pressure in pressures:
            dirname = (
                f"zeo_ff_{zeo_ff_name}_T_{temperature:0.1f}_P_{pressure:0.2f}".replace(
                    " ", "_"
                ).replace("/", "-")
            )
            thermo = ThermoProps(dirname + "/prod.out.prp")
            avg_n.append(np.mean(thermo.prop("Nmols_2")) / n_unitcells)
        avg_n = u.unyt_array(avg_n)
        ax.scatter(
            pressures.to_value("bar"),
            avg_n.value,
            marker="s",
            facecolors="none",
            edgecolors=prev_plt.get_facecolor()[0],
            label=f"{temperature:.0f} TraPPE",
            s=100,
        )

    ax.scatter(
        lit_298K[:, 0],
        lit_298K[:, 1],
        marker="1",
        label="Romanielo 2020 298 K",
        s=150,
    )
    ax.scatter(
        lit_309K[:, 0],
        lit_309K[:, 1],
        marker="2",
        label="Romanielo 2020 309 K",
        s=150,
    )

    ax.set_ylabel("$\mathregular{N_{methane}}$/unit cell", fontsize=16, labelpad=10)
    ax.set_xlabel("Pressure (bar)", fontsize=16, labelpad=10)
    ax.tick_params(which="both", direction="in", labelsize=14, top=True, right=True)
    ax.set_xscale("log")
    ax.set_ylim(0, 1.5)
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    fig.legend(fontsize=14, loc="upper left", bbox_to_anchor=(0.18, 0.95))
    fig.tight_layout(pad=2.0)

    fig.savefig(outdir + f"/isotherm.pdf")


if __name__ == "__main__":
    main()
