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

# Load in last frame of the GCMC simulation
gcmc_system = mb.load('last_frame.xyz')
n_water = (gcmc_system.n_particles - 2016) / 3

# Create separate mBuild molecules for water and graphene to save out gro or mol2 file
graphene = mb.Compound()
water = mb.Compound()
single_mol = mb.Compound()
gro_system = mb.Compound()

# Wrap water particles into residues
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

gro_system.add(mb.clone(graphene))
gro_system.add(mb.clone(water))

# Load graphene and water
graphene = recipes.GraphenePore(pore_width=1.5, pore_length=2.95, pore_depth=2.98, slit_pore_dim=2)
single_water = mb.load('structures/spce.mol2')

# Fix the box info
gcmc_system.periodicity = graphene.periodicity
gro_system.periodicity = graphene.periodicity
gro_system.save('coords.gro', combine='all', overwrite=True, residues=['C', 'SOL'])

# Create FF object and apply to graphene and water
ff = Forcefield('ffxml/C-spce.xml')
typed_graphene = ff.apply(graphene)
typed_water = ff.apply(single_water)

box_list = [gcmc_system]
species_list = [typed_graphene, typed_water]

# Specify mol to add in box
mols_to_add = [[0, 0]]
# Specify mols at start of the simulation
mols_in_boxes = [[1, int(n_water)]]
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

custom = {'coord_freq': 500}

# Run the NVT MC simulation
mc.run(
system=system,
moveset=moves,
run_type="equil",
run_length=1000000,
temperature=298.0 * u.K,
**custom
)

# Add coordinates into XYZ file
wrap_coords('nvt.out.xyz', graphene.periodicity*10)
