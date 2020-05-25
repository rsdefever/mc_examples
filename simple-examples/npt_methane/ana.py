import unyt as u
import numpy as np
import matplotlib.pyplot as plt

from mosdef_cassandra.analysis import ThermoProps

thermo = ThermoProps("npt.out.prp")

plt.plot(
    thermo.prop("MC_STEP").value,
    thermo.prop("Pressure").to_value(u.MPa),
)
plt.xlabel("MC Step", fontsize=24, labelpad=10)
plt.ylabel("Pressure (MPa)", fontsize=24, labelpad=10)
plt.xticks(np.arange(0,250001,75000), fontsize=16)
plt.yticks(fontsize=16)
plt.tight_layout()

plt.savefig("pressure.pdf")



