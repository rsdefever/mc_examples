import mdtraj as md
import os
from mtools.post_process import calc_msd
from mtools.gromacs.gromacs import unwrap_trj


def main():

    # Output
    pore_sizes_list = []
    n_ion_pair_list = []
    d_list = []

    for pore_size in pore_sizes:
        for n_ion_pair in n_ion_pairs:
            trr_path = (
                f"../md_pore/{pore_size.to_value('nm')}nm_{n_ion_pair}pairs/nvt.trr"
            )
            unwrap_path = f"../md_pore/{pore_size.to_value('nm')}nm_{n_ion_pair}pairs/nvt_unwrapped.trr"
            gro_path = (
                f"../md_pore/{pore_size.to_value('nm')}nm_{n_ion_pair}pairs/nvt.gro"
            )

            # If unwrap trajectory doesn't exist, call `unwrap_trj` function
            if not os.path.isfile(trr_path):
                unwrap_trj(trr_path)

            # Create MDTraj trajectory, chop off first 5000 frames
            trj = md.load(unwrap_path, top=gro_path)[5000:]
            water_trj = trj.atom_slice(trj.topology.select("water"))

            # Calculated 2-D MSD
            D, msd, x_fit, y_fit = calc_msd(water, dims=[1, 1, 0])

            # Append data to respective lists
            pore_sizes_list.append(pore_size)
            n_ion_pair_list.append(n_ion_pair)
            d_list.append(D)

    df = pd.DataFrame(columns=["pore_size_nm", "n_ion_pairs", "diffusivity_m^2_per_s"])
    df["pore_size_nm"] = np.array(pore_sizes_list)
    df["n_ion_pairs"] = np.array(n_ion_pair_list)
    df["diffusivity_m^2_per_s"] = np.array(d_list)
    df.to_csv("results_md_pore.csv")


if __name__ == "__main__":
    main()
