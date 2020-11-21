import unyt as u
import numpy as np
from mosdef_cassandra.analysis import ThermoProps
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, AutoMinorLocator
from matplotlib import rcParams
import seaborn

rcParams['font.sans-serif'] = "Arial"
rcParams['font.family'] = "sans-serif"

thermo_liq = ThermoProps("gemc.out.box1.prp")
thermo_vap = ThermoProps("gemc.out.box2.prp")

fig, ax = plt.subplots()
ax.plot(
    thermo_vap.prop("MC_STEP").value,
    thermo_vap.prop("Pressure").to_value(u.MPa),
    color="black",
)

ax.set_xlabel("MC Step", fontsize=18, labelpad=18)
ax.set_ylabel("Vapor pressure (MPa)", fontsize=18, labelpad=18)
ax.set_xticks(np.arange(0, 500001, 125000))
ax.tick_params(labelsize=16)
ax.tick_params(which="both", direction="in", right=True, top=True)
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.set_xlim(0, 500000)
ax.set_ylim(0, 5.0)
fig.tight_layout()

fig.savefig("pressure.pdf")

fig, ax = plt.subplots()
ax.plot(
    thermo_liq.prop("MC_STEP").value,
    thermo_liq.prop("Mass_Density").to_value('kg/m**3'),
    label="Liquid",
    color="#2aa134",
    linewidth=2.0
)
ax.plot(
    thermo_vap.prop("MC_STEP").value,
    thermo_vap.prop("Mass_Density").to_value('kg/m**3'),
    label="Vapor",
    color="#4b85b3",
    linewidth=2.0
)

ax.set_xlabel("MC Step", fontsize=18, labelpad=18)
ax.set_ylabel("Density (kg/m$\mathregular{^3}$)", fontsize=18, labelpad=18)
ax.set_xticks(np.arange(0, 500001, 125000))
ax.tick_params(labelsize=16)
ax.tick_params(which="both", direction="in", right=True, top=True)
ax.legend(loc="center right", fontsize=18)
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.set_xlim(0, 500000)
ax.set_ylim(0, 600)
fig.tight_layout()

fig.savefig("density.pdf")

fig, ax = plt.subplots()
ax.plot(
    thermo_liq.prop("MC_STEP").value,
    thermo_liq.prop("Nmols"),
    label="Liquid box",
    color="#2aa134",
    linewidth=2.5
)
ax.plot(
    thermo_vap.prop("MC_STEP").value,
    thermo_vap.prop("Nmols"),
    label="Vapor box",
    color="#4b85b3",
    linewidth=2.5
)

ax.set_xlabel("MC Step", fontsize=18, labelpad=12)
ax.set_ylabel("Number of molecules", fontsize=18, labelpad=5)
ax.set_xticks(np.arange(0, 500001, 125000))
ax.tick_params(labelsize=16)
ax.tick_params(which="both", direction="in", right=True, top=True)
ax.legend(loc="center left", fontsize=18)
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.set_ylim(0,1000)
fig.tight_layout()

fig.savefig("nmols.pdf")

fig, ax = plt.subplots()
ax.plot(
    thermo_liq.prop("MC_STEP").value,
    (thermo_liq.prop("Volume")**(1./3.)).to_value(u.angstrom),
    label="Liquid box",
    color="#2aa134",
    linewidth=2.5
)
ax.plot(
    thermo_vap.prop("MC_STEP").value,
    (thermo_vap.prop("Volume")**(1./3.)).to_value(u.angstrom),
    label="Vapor box",
    color="#4b85b3",
    linewidth=2.5
)

ax.set_xlabel("MC Step", fontsize=18, labelpad=18)
ax.set_ylabel("Box length (Angstrom)", fontsize=18, labelpad=18)
ax.set_xticks(np.arange(0, 500001, 125000))
ax.tick_params(labelsize=16)
ax.tick_params(which="both", direction="in", right=True, top=True)
ax.legend(loc="center left", fontsize=18)
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.set_ylim(0,200)
fig.tight_layout()
fig.savefig("boxl.pdf")

