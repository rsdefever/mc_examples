import mbuild as mb
import foyer
import mosdef_cassandra as mc
import numpy as np
import os
import unyt as u

from mbuild import recipes
from foyer import Forcefield
from mosdef_cassandra.utils.tempdir import temporary_cd
from cassandra_slitpore.utils.utils import translate_compound, wrap_coords
from utils import get_last_frame

# Create graphene system and atom type
graphene = recipes.GraphenePore(pore_width=1.5, pore_length=2.95, pore_depth=2.98, slit_pore_dim=2)
ff = Forcefield('ffxml/C-spce.xml')

typed_graphene = ff.apply(graphene)

# Create solvent
water = mb.load('structures/spce.mol2')
typed_water = ff.apply(water)

# Create box and species list
box_list = [graphene]
species_list = [typed_graphene, typed_water]

# Specify mol to add in box
mols_to_add = [[0, 296]]
# Specify mols at start of the simulation
mols_in_boxes = [[1, 0]]

# Create MC system
system = mc.System(box_list, species_list, mols_in_boxes=mols_in_boxes, mols_to_add=mols_to_add)
moves = mc.MoveSet("nvt", species_list)

# Set thermodynamic properties
thermo_props = [
    "energy_total",
    "pressure",
    "volume",
    "nmols",
    "mass_density",
]

mu = -40.0 * (u.kJ/u.mol)

custom_args = {
    "run_name": "equil",
    "charge_style": "ewald",
    "rcut_min": 0.2 * u.angstrom,
    "vdw_cutoff": 10.0 * u.angstrom,
    "units": "sweeps",
    "properties": thermo_props,
}


# Cassandra scales OK up to 4-8 cores
mc.run(
system=system,
moveset=moves,
run_type="equilibration",
run_length=100000,
temperature=298.0 * u.K,
**custom_args,
)

wrap_coords('equil.out.xyz', graphene.periodicity*10)
get_last_frame('equil.out.xyz', 100000)

# Load in last frame of the GCMC simulation
gcmc_system = mb.load('last_frame.xyz')

# Create separate mBuild molecules for water and graphene
graphene = mb.Compound()
water = mb.Compound()
single_mol = mb.Compound()
system = mb.Compound()

water_count = 0
for child in gcmc_system.children:
    if child.name == 'C':
        graphene.add(mb.clone(child))
    else:
        water_count += 1
        single_mol.add(mb.clone(child))
        if water_count % 3 == 0:
            single_mol.name = 'SOL'
            single_mol.generate_bonds('O', 'H', 0, 0.3)
            water.add(mb.clone(single_mol))
            single_mol = mb.Compound()

#water.add(mb.clone(graphene))
system.add(mb.clone(graphene))
system.add(mb.clone(water))
system.periodicity = np.array([2.9472, 2.97774175, 3.175])

# Create FF object and apply to graphene and water
ff = Forcefield('ffxml/C-spce.xml')

#typed_water = ff.apply(water)

system.save('coords.gro', combine='all', overwrite=True, residues=['C', 'SOL'])
