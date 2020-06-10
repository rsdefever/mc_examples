import mbuild as mb
import foyer
import mosdef_cassandra as mc
import numpy as np
import os
import unyt as u

from mbuild import recipes
from foyer import Forcefield
from mosdef_cassandra.utils.tempdir import temporary_cd
from cassandra_slitpore.utils.utils import translate_compound, wrap_coords, translate_back

# Load in last frame of the GCMC simulation
gcmc_system = mb.load('../last_frame.xyz')

# Create separate mBuild molecules for water and graphene to save out gro or mol2 file
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


# Load graphene and water
new_graphene = recipes.GraphenePore(pore_width=1.5, pore_length=2.95, pore_depth=2.98, slit_pore_dim=2)

# Fix the box info
water.periodicity = new_graphene.periodicity
graphene.periodicity = new_graphene.periodicity
translate_compound(water)
translate_compound(graphene)
translate_compound(gcmc_system)
translate_back(water)
translate_back(graphene)
translate_back(gcmc_system)

# Create FF object and apply to graphene and water
ff = Forcefield('../ffxml/C-spce.xml')
waterPM = ff.apply(water, residues='SOL')
graphenePM = ff.apply(graphene, residues='C')

systemPM = graphenePM + waterPM

systemPM.save('translate.gro', combine='all', overwrite=True)
systemPM.save('translate.mol2', overwrite=True)
systemPM.save('translate.top', combine='all', overwrite=True)
