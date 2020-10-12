import mbuild
import numpy as np
import os


from pkg_resources import resource_filename


def create_spce_water():
    """Generate a single water molecule with the SPC/E geometry
    Paper DOI: 10.1021/j100308a038
    Arguments
    ---------
    None
    Returns
    -------
    mbuild.Compound
    """
    OH_bl = 0.1  # nm
    HOH_angle = 109.47  # degrees
    water = mbuild.Compound(name="SOL")
    O = mbuild.Compound(name="O", pos=[0.0, 0.0, 0.0])
    H1 = mbuild.Compound(name="H", pos=[OH_bl, 0.0, 0.0])
    H2 = mbuild.Compound(
        name="H",
        pos=[
            OH_bl * np.cos(np.radians(HOH_angle)),
            OH_bl * np.sin(np.radians(HOH_angle)),
            0.0,
        ],
    )
    water.add([O, H1, H2])
    water.add_bond([O, H1])
    water.add_bond([O, H2])

    return water


# Instantiate a single water into the namespace
spce_water = create_spce_water()


def create_system(
    pore_width,
    n_ion_pairs,
    n_water,
    pve_ion = mbuild.Compound(name="Na"),
    nve_ion = mbuild.Compound(name="Cl"),
):
    """Create a system with a graphene pore filled with n_ion_pairs,
    (na, cl), and n_water spce water molecules.

    Parameters
    ----------
    pore_width : u.unyt_quantity (length)
        width of pore for restricted insertions
    n_ion_pairs : int
        number of na, cl ion pairs
    n_water : int
        number of water molecules
    pve_ion : mbuild.Compound, optional, default=Na
        positive ion
    nve_ion : mbuild.Compound, optional, default=Cl
        negative ion

    Returns
    -------
    filled_pore : mbuild.Compound
    """

    if n_ion_pairs == 0:
        filled_pore = mbuild.recipes.GraphenePoreSolvent(
            pore_length=2.3,
            pore_depth=2.2,
            pore_width=pore_width.to_value("nm"),
            n_sheets=3,
            slit_pore_dim=2,
            x_bulk=0,
            solvent=spce_water,
            n_solvent=n_water,
        )
    else:
        # Create pore system
        filled_pore = mbuild.recipes.GraphenePoreSolvent(
            pore_length=2.3,
            pore_depth=2.2,
            pore_width=pore_width.to_value("nm"),
            n_sheets=3,
            slit_pore_dim=2,
            x_bulk=0,
            solvent=[pve_ion, nve_ion, spce_water],
            n_solvent=[n_ion_pairs, n_ion_pairs, n_water],
        )

    # Translate to centered at 0,0,0 and make box larger in z
    box_center = filled_pore.periodicity/2.0
    filled_pore.translate(-box_center)

    return filled_pore


def load_final_xyz_frame(fname):
    """Return the final frame of a Cassandra .xyz file as an mbuild.Compound
    Assumes there is a .H file with the same name. E.g., if the .xyz file
    is 'equil.out.xyz', there should also be an 'equil.out.H' containing
    the box information.
    Arguments
    ---------
    fname : str
        path to of the xyz file
    """
    if not isinstance(fname, str):
        raise TypeError("'fname' must be a string")
    if fname[-4:] == ".xyz":
        fname = fname[:-4]

    data = []
    with open(fname + ".xyz") as f:
        for line in f:
            data.append(line.strip().split())

    for iline, line in enumerate(data):
        if len(line) > 0:
            if line[0] == "MC_STEP:":
                natom_line = iline - 1

    final_frame = data[natom_line + 2 :]
    natoms = int(data[natom_line][0])
    with open(fname + "-final.xyz", "w") as f:
        f.write(f"{natoms}\nAtoms\n")
        for coord in final_frame:
            f.write(
                "{}\t{}\t{}\t{}\n".format(
                    coord[0], coord[1], coord[2], coord[3],
                )
            )
    data = []
    with open(fname + ".H") as f:
        for line in f:
            data.append(line.strip().split())

    nspecies = int(data[-1][0])
    box_matrix = np.asarray(
        data[-(nspecies + 5) : -(nspecies + 2)], dtype=np.float32
    )
    assert box_matrix.shape == (3, 3)
    if np.count_nonzero(box_matrix - np.diag(np.diagonal(box_matrix))) > 0:
        raise ValueError("Only orthogonal boxes are currently supported")

    # If all is well load in the final frame
    frame = mbuild.load(fname + "-final.xyz")
    # mbuild.Compounds use nanometers!
    frame.periodicity = np.diagonal(box_matrix / 10.0)

    return frame


def get_ff(filename):
    """Get path to a file in ffxml directory
    """
    file_path = resource_filename('mc_examples',
            os.path.join('realistic_workflows/graphene_slitpore/resources/ffxml', filename)
    )

    return file_path
