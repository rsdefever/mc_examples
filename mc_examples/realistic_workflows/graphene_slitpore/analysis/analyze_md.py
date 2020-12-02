import mdtraj as md
import os
import unyt as u
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mtools.post_process import calc_msd
from mtools.gromacs.gromacs import unwrap_trj
from copy import deepcopy


def main():
    # Systems simulated
    pore_area = 2 * 22.104 * 21.270 * u.angstrom**2 # From .inp file
    pore_sizes = [1.0, 1.5, 2.0] * u.nm
    n_ion_pairs = [0, 4, 8]

    # Output
    pore_sizes_list = []
    n_ion_pair_list = []
    d_list = []

    fig, ax = plt.subplots()
    for pore_size in pore_sizes:
        for n_ion_pair in n_ion_pairs:
            print((pore_size, n_ion_pair))
            trr_path = (
                f"../md_pore/{pore_size.to_value('nm')}nm_{n_ion_pair}pairs/nvt.trr"
            )
            unwrap_path = f"../md_pore/{pore_size.to_value('nm')}nm_{n_ion_pair}pairs/nvt_unwrapped.trr"
            gro_path = (
                f"../md_pore/{pore_size.to_value('nm')}nm_{n_ion_pair}pairs/nvt.gro"
            )

            split_path = trr_path.split('/')
            unwrap_split = 'nvt_unwrapped.trr'
            unwrap_split_path = deepcopy(split_path)
            unwrap_split_path[-1] = unwrap_split
            if not os.path.isfile(unwrap_path):
                print("Unwrapping trajectory...")
                os.system('gmx trjconv -f {0} -o {1} -pbc nojump'.format(
                        '/'.join(split_path),
                        '/'.join(unwrap_split_path),
                ))

            # Create MDTraj trajectory, chop off first 5000 frames
            trj = md.load(unwrap_path, top=gro_path)[5000:]
            water_trj = trj.atom_slice(trj.topology.select("water"))

            # Calculated 2-D MSD
            print("Calculating MSD with mtools...")
            mtools_split_path = deepcopy(split_path)
            mtools_split_path[-1] = "mtools_msd.txt"
            #if not os.path.isfile("/".join(mtools_split_path)):
            D_bar, D_std, msd_bar = _run_multiple(water_trj, dims=[1, 1, 0])
            np.savetxt("/".join(mtools_split_path), np.transpose(np.vstack([trj.time[:1000], msd_bar])),
                       header=f"time(ps)\tmsd(nm^2)\t#D_bar={D_bar},D_std={D_std}")

            # Calculate MSD from gromacs
            print("Calculating MSD with GROMACS...")
            unwrap_split_path = deepcopy(split_path)
            unwrap_split_path[-1] = unwrap_split
            gro_split_path = deepcopy(split_path)
            gro_split_path[-1] = "nvt.gro"
            msd_split_path = deepcopy(split_path)
            for species, gmx_call in {"water": "water", "Na": 10, "Cl": 11}.items():
                msd_split_path[-1] = f"msd_{species}.xvg"
                if species in ("Na", "Cl") and n_ion_pair == 0:
                    continue
                os.system('echo {0} | gmx msd -f {1} -s {2} -o {3} -b 5000 --lateral z'.format(
                           gmx_call,
                           '/'.join(unwrap_split_path),
                           '/'.join(gro_split_path),
                           '/'.join(msd_split_path),
                ))
    
                # Append data to respective lists
                pore_sizes_list.append(pore_size)
                n_ion_pair_list.append(n_ion_pair)
                d_list.append(D_bar)
    
    df = pd.DataFrame(columns=["pore_size_nm", "n_ion_pairs", "diffusivity_m^2_per_s"])
    df["pore_size_nm"] = np.array(pore_sizes_list)
    df["n_ion_pairs"] = np.array(n_ion_pair_list)
    df["diffusivity_m^2_per_s"] = np.array(d_list)
    df.to_csv("results_md_pore.csv")

def _run_multiple(trj, dims):
    D_pop = list()
    msd_pop = list()
    for start_frame in np.linspace(0, trj.n_frames, num=200, dtype=np.int):
        end_frame = start_frame + 1000
        if end_frame < trj.n_frames:
            chunk = trj[start_frame:end_frame]
            print('\t\t\t...frame {} to {}'.format(start_frame, end_frame))
            try:
                D, msd, x_fit, y_fit = calc_msd(chunk, dims=dims)
                D_pop.append(D)
                msd_pop.append(msd)
            except TypeError:
                import pdb
                pdb.set_trace()
        else:
            continue
    D_bar = np.mean(D_pop)
    D_std = np.std(D_pop)
    msd_bar = np.mean(msd_pop, axis=0)
    return D_bar, D_std, msd_bar


if __name__ == "__main__":
    main()
