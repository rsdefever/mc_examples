import mbuild
import foyer
import unyt as u
import mosdef_cassandra as mc

# Load the force field from foyer
ff = foyer.forcefields.load_OPLSAA()

# Create an ethane molecule from the SMILES string
ethane = mbuild.load("CC", smiles=True)

# Apply the force field
ethane_ff = ff.apply(ethane)

# Initial vapor box based IG law @ 250 K, 10 bar
N_vapor = 50
pressure = 10.0 * u.bar
temperature = 250 * u.K
vapor_volume = N_vapor * u.kb * temperature / pressure
box_length = (vapor_volume ** (1./3.)).to_value(u.nm)
box = mbuild.Box([box_length, box_length, box_length])
vapor_box = mbuild.fill_box(ethane, n_compounds=N_vapor, box=box)

# Initial liquid box from existing configuration
liquid_box = mbuild.load("liqbox_equil.gro")
N_liquid = len([child for child in liquid_box.children])

# Create the System and MoveSet
box_list = [liquid_box, vapor_box]
species_list = [ethane_ff]
mols_in_boxes = [[N_liquid], [N_vapor]]

system = mc.System(
    box_list,
    species_list,
    mols_in_boxes=mols_in_boxes,
)
moveset = mc.MoveSet("gemc", species_list)

# Run the Monte Carlo simulation
mc.run(
    system=system,
    moveset=moveset,
    run_type="equilibration",
    run_length=500000,
    temperature=temperature,
    vdw_cutoff=10.0 * u.angstrom,
    charge_cutoff=10.0 * u.angstrom,
)


