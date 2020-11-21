import mbuild
import foyer
import unyt as u
import mosdef_cassandra as mc

# Load the force field from foyer
ff = foyer.forcefields.load_TRAPPE_UA()

# Create a single site methane molecule
methane = mbuild.Compound(name="_CH4")

# Apply the force field
methane_ff = ff.apply(methane)

# Create an empty simulation box
box = mbuild.Box([3.,3.,3.])

# Create the System and MoveSet
box_list = [box]
species_list = [methane_ff]
mols_to_add = [[100]]

system = mc.System(
    box_list, species_list, mols_to_add=mols_to_add
)
moveset = mc.MoveSet("npt", species_list)

# Run the Monte Carlo simulation
mc.run(
    system=system,
    moveset=moveset,
    run_type="equilibration",
    run_length=300000,
    temperature=240 * u.K,
    pressure=10.0 * u.bar,
)

