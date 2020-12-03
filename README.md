# mc_examples
This repository contains a collection of examples for the MoSDeF Cassandra Monte Carlo package. With the exception of the final slit pore example, all of the examples can be performed on a modern desktop or laptop in a reasonable period of time (< $\sim$1 hour).

## Helpful links
* MoSDeF Cassandra: [Docs](https://mosdef-cassandra.readthedocs.io), [GitHub](https://github.com/maginngroup/mosdef_cassandra)
* Cassandra: [Docs](https://cassandra-mc.readthedocs.io), [GitHub](https://github.com/maginngroup/cassandra)
* [MoSDeF tools](https://mosdef.org)
* [TRUE simulations](https://doi.org/10.1080/00268976.2020.1742938)
* [Installing Docker](https://docs.docker.com/get-docker/)
* [Installing Miniconda](https://docs.conda.io/en/latest/miniconda.html)

## Installation instructions

There are two options for installing the packages required to run the simulations in this repository: (1) via conda, or (2) via docker. If you are working in an HPC environment, you may need to use [singularity](https://sylabs.io/docs/#singularity) rather than docker. Installation instructions for miniconda and docker can be found in the "helpful links" section above.

### conda

Make sure you have conda installed. Then, clone the repository to your machine:

```
git clone https://github.com/rsdefever/mc_examples.git
```

Create a new conda environment (``mc-ex``) and install the required packages:

```
conda create --name mc-ex --file mc_examples/requirements.txt
```

Finally, install the ``mc_examples`` via pip:

```
cd mc_examples/
pip install .
```

This final step is only strictly necessary if you plan on running the slit pore example.

### docker

Make sure you have docker installed. Get the docker image from dockerhub with the following:

```
docker pull rsdefever/mc_examples:latest
```

## Running the simple examples

### With conda-based installation
You may run the examples directly within the repository or in another location on your machine. If the ``mc_examples`` is installed at ``PATH_TO_MC_EXAMPLES``, and you want to run the NpT example under ``PATH_TO_WORKSPACE/npt``, you would do the following

```
conda activate mc-ex
cd PATH_TO_WORKSPACE
mkdir npt
cd npt/
python PATH_TO_MC_EXAMPLES/mc_examples/simple_examples/npt_methane/npt.py
python PATH_TO_MC_EXAMPLES/mc_examples/simple_examples/npt_methane/analyze.py
```

To run the GEMC example:

```
conda activate mc-ex
cd PATH_TO_WORKSPACE
mkdir gemc
cd gemc/
python PATH_TO_MC_EXAMPLES/mc_examples/simple_examples/gemc_ethane/gemc.py
python PATH_TO_MC_EXAMPLES/mc_examples/simple_examples/gemc_ethane/analyze.py
```

### With docker-based installation

First, change to the directory on your local machine where you want to store the simulation output. This is notated below as ``PATH_TO_WORKSPACE``. Then, we call the appropriate script from within the docker container.

To run the NpT example:

```
cd PATH_TO_WORKSPACE/
mdkir npt
cd npt/
docker run -t --mount type=bind,source=$(pwd),target=/workspace rsdefever/mc_examples:latest "python /software/mc_examples/mc_examples/simple_examples/npt_methane/npt.py"
docker run -t --mount type=bind,source=$(pwd),target=/workspace rsdefever/mc_examples:latest "python /software/mc_examples/mc_examples/simple_examples/npt_methane/analyze.py"
```

To run the NpT example:

```
cd PATH_TO_WORKSPACE/
mdkir gemc
cd gemc/
docker run -t --mount type=bind,source=$(pwd),target=/workspace rsdefever/mc_examples:latest "python /software/mc_examples/mc_examples/simple_examples/gemc_ethane/npt.py"
docker run -t --mount type=bind,source=$(pwd),target=/workspace rsdefever/mc_examples:latest "python /software/mc_examples/mc_examples/simple_examples/gemc_ethane/analyze.py"
```


## Running the zeolite adsorption example

Once again, we assume ``mc_examples`` is installed at ``PATH_TO_MC_EXAMPLES``, and you want to run the examples under ``PATH_TO_WORKSPACE/``.

### With conda-based installation

```
conda activate mc-ex
cd PATH_TO_WORKSPACE
mkdir zeolite
cd zeolite
cp -r PATH_TO_MC_EXAMPLES/mc_examples/realistic_workflows/zeolite_adsorption/resources .
python PATH_TO_MC_EXAMPLES/mc_examples/realistic_workflows/zeolite_adsorption/run_fluid.py
python PATH_TO_MC_EXAMPLES/mc_examples/realistic_workflows/zeolite_adsorption/analyze_fluid.py
python PATH_TO_MC_EXAMPLES/mc_examples/realistic_workflows/zeolite_adsorption/run_adsorption.py
python PATH_TO_MC_EXAMPLES/mc_examples/realistic_workflows/zeolite_adsorption/analyze_adsorption.py
```

The output from the analysis scripts will be located in a directory named ``plots/``.

### With docker-based installation

```
cd PATH_TO_WORKSPACE/
mdkir zeolite
cd zeolite/

```

## Running the slit pore example

### With conda-based installation

### With docker-based installation