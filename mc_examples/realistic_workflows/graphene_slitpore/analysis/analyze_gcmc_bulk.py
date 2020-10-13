import unyt as u
import numpy as np
import pandas as pd

from mosdef_cassandra.analysis import ThermoProps


def main():

    # Input conditions
    temperature = 298.0 * u.K
    mus = np.arange(-58, -42, 2) * u.Unit("kJ/mol")

    # Output
    pressures = []

    for mu in mus:
        dirname = f"T_{temperature:0.1f}_mu_{mu:.1f}".replace(
            " ", "_"
        ).replace(
            "/", "-"
        )
        thermo_path = "../gcmc_bulk/" + dirname + "/prod.out.prp"
        thermo = ThermoProps(thermo_path)
        pressures.append(thermo.prop("Pressure").mean())

    pressures = u.unyt_array(pressures)

    df = pd.DataFrame(
        columns=["mu-cassandra_kJmol", "pressure_bar"]
    )
    df["mu-cassandra_kJmol"] = mus.to_value("kJ/mol")
    df["pressure_bar"] = pressures.to_value("bar")
    df.to_csv("results_gcmc_bulk.csv")


if __name__ == "__main__":
    main()
