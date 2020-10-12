import mbuild
import foyer
import unyt as u
import mosdef_cassandra as mc

from mc_examples.realistic_workflows.graphene_slitpore.utils import spce_water, load_final_xyz_frame, get_ff


def run_gcmc(
    filled_pore,
    pore_width,
    temperature,
    mu,
    nsteps_nvt,
    nsteps_gcmc,
    pve_ion = mbuild.Compound(name="Na"),
    nve_ion = mbuild.Compound(name="Cl"),
    **custom_args,
):
    """Run desorption simulation at the specified temperature
    and chemical potential
    Parameters
    ----------
    filled_pore : porebuilder.GraphenePoreSolvent
        pore filled with water and (optional) ions
    pore_width : u.unyt_quantity (length)
        width of pore for restricted insertions
    temperature: u.unyt_quantity (temperature)
        desired temperature
    mu : u.unyt_quantity (energy)
        desired chemical potential for water
    nsteps_nvt : int
        number of MC steps for NVT equilibration
    nsteps_gcmc : int
        number of MC steps for GCMC simulation
    pve_ion : mbuild.Compound, optional, default=Na
        positive ion
    nve_ion : mbuild.Compound, optional, default=Cl
        negative ion
    custom_args : opt, additional keyword arguments
        provide additional custom keyword arguments to MoSDeF Cassandra
        and override the default values
    Returns
    -------
    None: runs simulation
    """

    # Load foyer ff
    ff = foyer.Forcefield(get_ff("pore-spce-jc.xml"))

    # Extract just the pore and apply ff
    empty_pore = filled_pore.children[0]
    typed_pore = ff.apply(empty_pore)

    # Create a water molecule with the spce geometry and apply ff
    typed_water = ff.apply(spce_water)
    typed_pve = ff.apply(pve_ion)
    typed_nve = ff.apply(nve_ion)

    # Determine the number of waters in the pore
    n_water = len([child for child in filled_pore.children if child.name=="SOL"])
    n_ion_pairs = int((len([child for child in filled_pore.children])-1-n_water)/2)

    # Create box and species list
    box_list = [filled_pore]
    species_list = [typed_pore, typed_pve, typed_nve, typed_water]

    # Specify mols at start of the simulation
    mols_in_boxes = [[1, n_ion_pairs, n_ion_pairs, n_water]]

    # Create MC system
    system = mc.System(box_list, species_list, mols_in_boxes=mols_in_boxes)
    moves = mc.MoveSet("nvt", species_list)

    # Set move probabilities
    moves.prob_translate = 0.5
    moves.prob_rotate = 0.5
    moves.prob_regrow = 0.0

    # Set thermodynamic properties
    thermo_props = [
        "energy_total",
        "energy_intervdw",
        "energy_interq",
        "nmols",
    ]

    default_args = {
        "run_name" : "nvt",
        "cutoff_style": "cut",
        "charge_style": "ewald",
        "rcut_min": 0.5 * u.angstrom,
        "vdw_cutoff": 9.0 * u.angstrom,
        "charge_cutoff": 9.0 * u.angstrom,
        "properties": thermo_props,
        "angle_style": ["harmonic", "harmonic", "harmonic", "fixed"],
        "coord_freq": 100000,
        "prop_freq": 1000,
    }

    custom_args = { **default_args, **custom_args}

    # Run NVT equilibration
    mc.run(
        system=system,
        moveset=moves,
        run_type="equilibration",
        run_length=nsteps_nvt,
        temperature=temperature,
        **custom_args,
    )

    # Create MC system
    equilibrated_box = load_final_xyz_frame("nvt.out.xyz")
    box_list = [equilibrated_box]
    system = mc.System(box_list, species_list, mols_in_boxes=mols_in_boxes)
    moves = mc.MoveSet("gcmc", species_list)

    # Set move probabilities
    moves.prob_translate = 0.25
    moves.prob_rotate = 0.25
    moves.prob_insert = 0.25
    moves.prob_regrow = 0.0

    # Specify the restricted insertion
    restricted_type = [[None, None, None, "slitpore"]]
    restricted_value = [[None, None, None, 0.5 * pore_width]]
    moves.add_restricted_insertions(
        species_list, restricted_type, restricted_value
    )

    # Run GCMC
    custom_args["run_name"] = "gcmc"
    mc.run(
        system=system,
        moveset=moves,
        run_type="equilibration",
        run_length=nsteps_gcmc,
        temperature=temperature,
        chemical_potentials=["none", mu],
        **custom_args,
    )

