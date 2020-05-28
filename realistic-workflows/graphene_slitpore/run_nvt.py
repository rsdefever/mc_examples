import mbuild as mb
import foyer
import mosdef_cassandra as mc
import numpy as np
import os
import unyt as u

from mbuild import recipes
from foyer import Forcefield
from mosdef_cassandra.utils.tempdir import temporary_cd
from cassandra_slitpore.utils.utils import translate_compound

# Load in last frame of the GCMC simulation
gcmc_system = mb.load('last_frame.xyz')

# Create separate mBuild molecules for water and graphene
graphene = mb.Compound()
water = mb.Compound()
single_mol = mb.Compound()

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

water.add(mb.clone(graphene))
water.periodicity = np.array([2.9472, 2.97774175, 23.175])
# Create solvent
single_water = mb.load('structures/spce.mol2')

# Create FF object and apply to graphene and water
ff = Forcefield('ffxml/C-spce.xml')
typed_graphene = ff.apply(graphene)
typed_water = ff.apply(single_water)
#typed_water = ff.apply(water.to_parmed(infer_residues=True))

# Create box and species list
#box_list = [graphene, water]
box_list = [water]
#species_list = [typed_graphene, typed_water]
#species_list = [i for i in typed_water.residues]
#species_list.insert(0, typed_graphene)
species_list = [typed_graphene, typed_water]

# Specify mol to add in box
mols_to_add = [[0, 0]]
#mols_to_add = [0] * (len(typed_water.residues)+1)
#mols_to_add = [mols_to_add]
# Specify mols at start of the simulation
mols_in_boxes = [[1,297]]
#mols_in_boxes = [1] * (len(typed_water.residues)+1)
#mols_in_boxes = [mols_in_boxes]
import pdb; pdb.set_trace()

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

# Cassandra scales OK up to 4-8 cores
mc.run(
system=system,
moveset=moves,
run_type="production",
run_length=2000,
temperature=298.0 * u.K,
)
