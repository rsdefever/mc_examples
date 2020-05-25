import mbuild
import foyer
import unyt as u
import mosdef_cassandra as mc

# Load the force field from foyer
ff = foyer.forcefields.load_TRAPPE_UA()

# Create a single site methane
methane = mbuild.Compound(name="_CH4")

# Apply the force field
methane_ff = ff.apply(methane)

# Initial vapor box based IG law @ 110 K, 1 bar
N_vapor = 200
pressure = 1.0 * u.bar
temperature = 110 * u.K
vapor_volume = N_vapor * u.kb * temperature / pressure
box_length = (vapor_volume ** (1./3.)).to_value(u.nm)
vapor_box = mbuild.Box([box_length, box_length, box_length])

# Initial liquid box from existing configuration
liquid_box = mbuild.load("liqbox_equil.gro")
N_liquid = liquid_box.n_particles

# Create the System and MoveSet
box_list = [liquid_box, vapor_box]
species_list = [methane_ff]
mols_in_boxes = [[N_liquid], [0]]
mols_to_add = [[0], [N_vapor]]

system = mc.System(
    box_list,
    species_list,
    mols_in_boxes=mols_in_boxes,
    mols_to_add=mols_to_add
)
moveset = mc.MoveSet("gemc", species_list)

# Run the Monte Carlo simulation
mc.run(
    system=system,
    moveset=moveset,
    run_type="equilibration",
    run_length=500000,
    temperature=temperature,
    charge_style="none",
    vdw_cutoff=14.0 * u.angstrom,
)

