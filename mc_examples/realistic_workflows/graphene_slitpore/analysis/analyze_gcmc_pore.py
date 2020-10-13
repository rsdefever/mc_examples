import unyt as u
import numpy as np
import pandas as pd

from mosdef_cassandra.analysis import ThermoProps


def main():

    # Systems simulated
    pore_area = 2 * 22.104 * 21.270 * u.angstrom**2 # From .inp file
    pore_sizes = [1.0, 1.5, 2.0] * u.nm
    n_ion_pairs = [0, 4, 8]

    # Output
    nmols_list = []
    pore_sizes_list = []
    n_ion_pair_list = []

    for pore_size in pore_sizes:
        for n_ion_pair in n_ion_pairs:
            thermo_path = f"../gcmc_pore/{pore_size.to_value('nm')}nm_{n_ion_pair}pairs/gcmc.out.prp"
            thermo = ThermoProps(thermo_path)
            nmols_list.append(thermo.prop("Nmols_4", start=20000000).mean())
            pore_sizes_list.append(pore_size)
            n_ion_pair_list.append(n_ion_pair)

    df = pd.DataFrame(
        columns=["pore_size_nm", "n_ion_pairs", "nmols", "nmols_per_nm^2"]
    )
    df["pore_size_nm"] = np.array(pore_sizes_list)
    df["n_ion_pairs"] = np.array(n_ion_pair_list)
    df["nmols"] = np.array(nmols_list)
    df["nmols_per_nm^2"] = np.array(nmols_list) / pore_area.to_value(u.nm**2)
    df.to_csv("results_gcmc_pore.csv")


if __name__ == "__main__":
    main()
