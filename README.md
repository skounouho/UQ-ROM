# Reduced-Order Modeling for LAMMPS Using Proper Orthogonal Decomposition

This collection of fixes allow for the creation of reduced order models in LAMMPS using POD linear subspace methods.

## Installation

Download the files in `src` and insert them in the `src` folder in your LAMMPS installation. You will need to install [Eigen](https://eigen.tuxfamily.org) to use `fix rob`. See the [LAMMPS documentation](https://eigen.tuxfamily.org) for how to compile LAMMPS with Eigen.

## Background

The goal of proper orthogonal decomposition is to take a complex system of seemingly random vectors and extract some kind of order from the chaos. This is done by modeling the trajectory of each particle as a function of a spatially-dependent function and a time-dependent coefficient.

$$
\pmb{q}(\pmb{x},t) = \sum^{\infty}_{k=1} \pmb{\Phi}_k(\pmb{x}) a_k(t)
$$

The reduced order model works by taking $N_s$ atom displacement snapshots of a full atomistic simulation and performing POD using singular value decomposition. The first matrix $[\Phi]$ in the resulting decomposition can be used as a linear approximation of the spatially-dependent functions, while the time-dependent coefficients become the reduced variable $\pmb{y}$.

We can create a reduced order system by projection from the physical space $\pmb{q}$ to the reduced space.

$$
\pmb{y}(t) = [\Phi]^T (\pmb{q}(t) - \pmb{q}(0))
$$

For a longer primer on the POD method, see [Weiss 2019](https://doi.org/10.2514/6.2019-3333). For a full mathematical proof of the method, see [Gubisch and Volkwein 2017](https://doi.org/10.1137/1.9781611974829.ch1).

## Usage

### fix nve/rom command

```
fix ID group-ID nve/rom order robfile
```
* ID, group-ID are documented in [fix](https://docs.lammps.org/fix.html) command
* nve/rom = style name of this fix command
* order = model order, which corresponds to the number of columns to be read from the reduced-order basis
* robfile = filepath to the input file containing the reduced-order basis

### fix nvt/rom, npt/rom, nph/rom commands
```
fix ID group-ID style_name keyword value ...
```
* ID, group-ID are documented in [fix](https://docs.lammps.org/fix.html) command
* style_name = *nvt/rom*, *npt/rom* or *nph/rom*
* one keyword must be appended
```
keyword = model
values = order robfile
  order = model order, which corresponds to the number of columns to be read from the
      reduced-order basis
  robfile = filepath to the input file containing the reduced-order basis
```

The `rom` fixes run a reduced-order simulation in the chosen ensemble. The code reads the reduced-order basis, converts between physical and reduced quantities and calculate time integration in the reduced-order space. The reduced model *should* be compatible with running on multiple processors using MPI.

NOTE: To use the `rom` command with the NVT, NPT or NPH ensembles, make sure to replace the `fix_nh_rom.cpp` file in the LAMMPS `src` directory. The file includes a few lines of code for parsing the keyword 'model'.

### fix rob command

```
fix ID group-ID rob N order robfile
```
* ID, group-ID are documented in [fix](https://docs.lammps.org/fix.html) command
* rob = style name of this fix command
* N = takes snapshots every N timesteps
* order = model order, which corresponds to the number of columns to be written from the reduced-order basis
* robfile = filepath to the output file for the reduced-order basis

The `rob` fix generates a reduced order basis using the POD method. The code assembles an array of snapshots during the run, computes the reduced-order basis and writes the results to a user-designated file. The output file format is space-delimited, with the singular values listed in decreasing order at the bottom of the file.

This command can be used with the LAMMPS [rerun](https://docs.lammps.org/rerun.html) command, as long as the timestep is set to a multiple of the timesteps in the dump files. **To decrease the memory usage during the initial run of the full-simulation, generating the reduced-order basis from a rerun is recommended.** Generating from a rerun also allows the user to create an basis for a simulation run on multiple processors.
