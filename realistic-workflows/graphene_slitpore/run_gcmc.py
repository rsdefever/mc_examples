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

# Translate the graphene structure
translate_compound(graphene)

typed_graphene = ff.apply(graphene)

# Create solvent
water = mb.load('structures/spce.mol2')
typed_water = ff.apply(water)

# Create box and species list
box_list = [graphene]
species_list = [typed_graphene, typed_water]

# Specify mol to add in box
mols_to_add = [[0, 10]]
# Specify mols at start of the simulation
mols_in_boxes = [[1, 0]]

# Create MC system
system = mc.System(box_list, species_list, mols_in_boxes=mols_in_boxes, mols_to_add=mols_to_add)
moves = mc.MoveSet("gcmc", species_list)

# Set move probabilities
moves.prob_translate = 0.35
moves.prob_rotate = 0.20

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

# Specify the restricted insertion
restricted_type = [[None, 'slitpore']]
restricted_value = [[None, 7.5 * u.angstrom]]
moves.add_restricted_insertions(species_list,
                               restricted_type,
                               restricted_value)

mc.run(
system=system,
moveset=moves,
run_type="equilibration",
run_length=300000,
temperature=298.0 * u.K,
chemical_potentials=["none", mu],
**custom_args,
)

# Write out trajectory with box coordinates
#wrap_coords('equil.out.xyz', graphene.periodicity*10)
# Get last frame of xyz trajectory
get_last_frame('equil.out.xyz', 300000)
wrap_coords('last_frame.xyz', graphene.periodicity*10)
