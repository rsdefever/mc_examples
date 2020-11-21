import mbuild
import foyer

# Load the force field from foyer
ff = foyer.forcefields.load_OPLSAA()

ethane = mbuild.load("CC", smiles=True)

# Create an empty simulation box
box = mbuild.fill_box(ethane, n_compounds=450, density=500)

# Apply the force field
typed = ff.apply(box)

typed.save("in.gro", overwrite=True)
typed.save("topol.top", overwrite=True)
