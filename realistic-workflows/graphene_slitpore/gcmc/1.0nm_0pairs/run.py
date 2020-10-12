import unyt as u
import sys


sys.path.append("../../resources")
sys.path.append("../")
from utils import create_system
from runners import run_gcmc


def main():

    pore_width = 1.0 * u.nm
    n_ion_pairs = 4
    n_water = 75

    temperature = 298 * u.K
    mu = -48 * u.kJ / u.mol
    nsteps_nvt = 5000000
    nsteps_gcmc = 45000000

    filled_pore = create_system(pore_width, n_ion_pairs, n_water)

    # Run the MC simulations
    run_gcmc(
        filled_pore,
        pore_width,
        temperature,
        mu,
        nsteps_nvt,
        nsteps_gcmc,
        seeds = [3005929952, 14958885759]
    )



if __name__ == "__main__":
    main()
