# mc_examples
This repository contains a collection of examples for the MoSDeF Cassandra Monte Carlo package. With the exception of the final slit pore example, all of the examples can be performed on a modern desktop or laptop in a reasonable period of time (< 1-2 hours).

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
conda create --name mc-ex --file mc_examples/requirements.txt -c conda-forge -c bioconda
```

Finally, install the ``mc_examples`` via pip:

```
cd mc_examples/
pip install .
```

Note if you are running OSX, the latest available version of ``gromacs`` from ``bioconda`` is ``2019.1``. The simulations reported in the manuscript were performed with version ``2020``.

### docker

Make sure you have docker installed. You may need to increase the memory available to containers. Allowing the container access to at least 8 GB of RAM is recommended. Get the docker image from dockerhub with the following:

```
docker pull rsdefever/mc_examples:latest
```

## Running the simple examples

### With conda-based installation
You may run the examples directly within the repository or in another location on your machine. If the ``mc_examples`` is installed at ``PATH_TO_MC_EXAMPLES``, and you want to run the NpT example under ``PATH_TO_WORKSPACE/npt``, you would do the following

```
conda activate mc-ex
export OMP_NUM_THREADS=4
cd PATH_TO_WORKSPACE
mkdir npt
cd npt/
python PATH_TO_MC_EXAMPLES/mc_examples/simple_examples/npt_methane/npt.py
python PATH_TO_MC_EXAMPLES/mc_examples/simple_examples/npt_methane/analyze.py
```

Here, the ``export OMP_NUM_THREADS=4`` line sets Cassandra to use 4 Open MP threads. You can select a different value if you desire. Values of 1, 2, 4, or 8 are recommended.

To run the GEMC example:

```
conda activate mc-ex
export OMP_NUM_THREADS=4
cd PATH_TO_WORKSPACE
mkdir gemc
cd gemc/
cp PATH_TO_MC_EXAMPLES/mc_examples/simple_examples/gemc_ethane/liqbox_equil.gro .
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

To run the GEMC example:

```
cd PATH_TO_WORKSPACE/
mdkir gemc
cd gemc/
cp PATH_TO_MC_EXAMPLES/mc_examples/simple_examples/gemc_ethane/liqbox_equil.gro .
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
docker run -t --mount type=bind,source=$(pwd),target=/workspace rsdefever/mc_examples:latest "python /software/mc_examples/mc_examples/realistic_workflows/zeolite_adsorption/run_fluid.py"
docker run -t --mount type=bind,source=$(pwd),target=/workspace rsdefever/mc_examples:latest "python /software/mc_examples/mc_examples/realistic_workflows/zeolite_adsorption/analyze_fluid.py"
docker run -t --mount type=bind,source=$(pwd),target=/workspace rsdefever/mc_examples:latest "python /software/mc_examples/mc_examples/realistic_workflows/zeolite_adsorption/run_adsorption.py"
docker run -t --mount type=bind,source=$(pwd),target=/workspace rsdefever/mc_examples:latest "python /software/mc_examples/mc_examples/realistic_workflows/zeolite_adsorption/analyze_adsorption.py"
```

## Running the slit pore example

The slit pore example contains longer simulations that likely require submitting to a job scheduler. This will mean the exact instructions are cluster-dependent. Therefore, we provide somewhat less specific instructions and focus on the `conda` installation. In principle, these simulations can be run from within the docker container as well.

Inside the `graphene_slitpore` directory, you will find three subdirectories for running simulations: `gcmc_bulk`, `gcmc_pore`, and `md_pore`. The `gcmc_bulk` contains `run.py`, which will perform GCMC simulations of SPC/E water vapor at a range of chemical potentials. Here is an example of a job submission script that would activate the `conda` environment, set the `OMP_NUM_THREADS` variable to `4`, and run the script (for a UGE scheduler):

```
#!/bin/bash
#$ -N water
#$ -pe smp 4
#$ -r n
#$ -m ae

conda activate mc-ex
export OMP_NUM_THREADS=4

python run.py
```

Inside the `analysis` directory, you will find the `analyze_gcmc_bulk.py` script, which you can run after the `gcmc_bulk` simulations are complete.

The GCMC simulations water and ions in the graphene pores are found under the `gcmc_pore` subdirectory. Each pore/number of ion pair combination is contained within a different subdirectory. Each contains a `run.py`, which once again can be submitted to a job scheduler. Once the simulations have completed, the `analyze_gcmc_pore.py` script can be used to perform the analysis.

Finally, the `md_pore` subdirectory contains the code required to run the MD simulations. This will require a `gromacs` installation with an executable named `gmx` and accessible on your `PATH`. If your `gromacs` executable has a different name, you can edit the `_gromacs_str` function under `md_pore/runners.py`.  Similar to the structure of the `gcmc_pore` subdirectory, each pore/number of ion pair combinations is contained within a different subdirectory in `md_pore`.  Each subdirectory contains a `run.py` script which can similarly be submitted to a job scheduler.  Once the simulations have completed, analysis of the lateral mean squared displacements (MSDs) through the GROMACS built-in tools can be computed with the `analyze_md.py` script contained in the `analysis` directory.  After the MSDs have been computed, the self-diffusion coefficients can be calculated with `fit.py` and visualized with `plot.py` in the `analysis` subdirectory contained within `md_pore`.  Upon running `fit.py`, the self-diffusion coefficients are also written to a `csv` file for future reference.
