import unyt as u
import matplotlib.pyplot as plt
import numpy as np

from mosdef_cassandra.analysis import ThermoProps

thermo_liq = ThermoProps("gemc.out.box1.prp")
thermo_vap = ThermoProps("gemc.out.box2.prp")

plt.plot(
    thermo_vap.prop("MC_STEP").value,
    thermo_vap.prop("Pressure").to_value(u.kPa),
    color="#2aa134",
)

plt.xlabel("MC Step", fontsize=24, labelpad=20)
plt.ylabel("Vapor pressure (kPa)", fontsize=24, labelpad=20)
plt.xticks(np.arange(0, 500001, 125000), fontsize=16)
plt.yticks(fontsize=16)
plt.tight_layout()

plt.savefig("pressure.pdf")
plt.close()

plt.plot(
    thermo_liq.prop("MC_STEP").value,
    thermo_liq.prop("Mass_Density").to_value('kg/m**3'),
    label="Liquid",
    color="#2aa134",
)
plt.plot(
    thermo_vap.prop("MC_STEP").value,
    thermo_vap.prop("Mass_Density").to_value('kg/m**3'),
    label="Vapor",
    color="#4b85b3"
)

plt.xlabel("MC Step", fontsize=24, labelpad=20)
plt.ylabel("Density (kg/m$\mathregular{^3}$)", fontsize=24, labelpad=20)
plt.xticks(np.arange(0, 500001, 125000), fontsize=16)
plt.yticks(fontsize=16)
plt.legend(loc="center right", fontsize=24)
plt.tight_layout()

plt.savefig("density.pdf")
