import unyt as u
import matplotlib.pyplot as plt
from mosdef_cassandra.analysis import ThermoProps
from matplotlib.ticker import MultipleLocator
from matplotlib import rcParams

rcParams['font.sans-serif'] = "Arial"
rcParams['font.family'] = "sans-serif"


thermo = ThermoProps("npt.out.prp")

fig, ax = plt.subplots()
ax.plot(
    thermo.prop("MC_STEP").value,
    thermo.prop("Pressure").to_value(u.MPa),
    color="black",
)


ax.set_xlabel("MC Step", fontsize=18, labelpad=15)
ax.set_ylabel("Pressure (MPa)", fontsize=18, labelpad=18)


ax.set_xlim(0,300000)
ax.set_ylim(0,14)
ax.xaxis.set_major_locator(MultipleLocator(75000))
ax.xaxis.set_minor_locator(MultipleLocator(25000))
ax.yaxis.set_major_locator(MultipleLocator(2))
ax.yaxis.set_minor_locator(MultipleLocator(1))
ax.tick_params(labelsize=16)
ax.tick_params(which="both", direction="in", right=True, top=True)
fig.tight_layout()

fig.savefig("pressure.pdf")
