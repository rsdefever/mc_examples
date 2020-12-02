import mbuild
import foyer
import os
import unyt as u
from mc_examples.realistic_workflows.graphene_slitpore.utils import (
    spce_water,
    create_system,
    get_ff,
)


def run_gromacs(
    filled_pore,
    pve_ion = mbuild.Compound(name="Na"),
    nve_ion = mbuild.Compound(name="Cl"),
):
    """Run MD simulation in GROMACS at the specified temperature
    Parameters
    ----------
    filled_pore : porebuilder.GraphenePoreSolvent
         pore filled with water and (optional) ions
    pve_ion : mbuild.Compound, optional, default=Na
        positive ion
    nve_ion : mbuild.Compound, optional, default=Cl
        negative ion
    Returns
    -------
    None: runs simulation
    """

    if not os.path.isfile("slitpore.gro"):
        # Create mb.Compound of water and ions in slitpore

        # Load foyer ff
        ff = foyer.Forcefield(get_ff("pore-spce-jc.xml"))

        # Extract just the pore and apply ff
        empty_pore = filled_pore.children[0]
        typed_pore = ff.apply(empty_pore)

        # Extract the electrolyte molecules
        electrolyte = filled_pore.children[1:]

        # Extract water and apply ff
        water = mbuild.Compound(
            [mbuild.clone(child) for child in electrolyte if child.name == "SOL"]
        )
        typed_electrolyte = ff.apply(water, residues="SOL")

        # Extract ions if there are any
        cations = mbuild.Compound(
            [
                mbuild.clone(child)
                for child in electrolyte
                if child.name == pve_ion.name
            ]
        )
        anions = mbuild.Compound(
            [
                mbuild.clone(child)
                for child in electrolyte
                if child.name == nve_ion.name
            ]
        )
        if len(cations.children) > 0:
            typed_cations = ff.apply(cations, residues=pve_ion.name)
            typed_anions = ff.apply(anions, residues=nve_ion.name)
            typed_electrolyte += typed_cations + typed_anions

        # Combine ParmEd structures
        typed_slitpore = typed_pore + typed_electrolyte

        # Save ParmEd structure to GROMACS files
        typed_slitpore.save("slitpore.gro", combine="all", overwrite=True)
        typed_slitpore.save("slitpore.top", combine="all", overwrite=True)

        # Create idx file
        cmd = f'printf "! r RES\nq\n" | gmx make_ndx -f slitpore.gro -o index.ndx'
        os.system(cmd)

    # Run energy minimization
    _gromacs_str("em", "slitpore")

    # Run NVT simulation
    _gromacs_str("nvt", "em")


def _gromacs_str(op_name, gro_name):
    mdp = f'../mdp_files/{op_name}'
    cmd = f'gmx grompp -f {mdp}.mdp -c {gro_name}.gro -p slitpore.top -n index.ndx -o {op_name}.tpr --maxwarn 1 && gmx mdrun -deffnm {op_name} -ntmpi 1'

    return os.system(cmd)
