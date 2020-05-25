import unyt as u
import matplotlib.pyplot as plt

from mosdef_cassandra.analysis import ThermoProps

thermo = ThermoProps("npt.out.prp")

plt.plot(
    thermo.prop("MC_STEP").value,
    thermo.prop("Pressure").to_value(u.MPa),
)
plt.xlabel("MC Step", fontsize=16, labelpad=10)
plt.ylabel("Pressure (MPa)", fontsize=16, labelpad=10)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.tight_layout()

plt.savefig("pressure.pdf")



