# Reduced-Order Modeling and Uncertainty Quantification for LAMMPS Using Proper Orthogonal Decomposition

This collection of fixes allow for the creation of reduced order models in LAMMPS using POD linear subspace methods.

## Installation

Download the files in `src` and insert them in the `src` folder in your LAMMPS installation. You will need to install [Eigen](https://eigen.tuxfamily.org) to use `fix rob`. I do not have a CMake file to do this, so you will need to symlink `Eigen` and `unsupported/Eigen` to the `src` directory in LAMMPS. Neither library has to be built separately.

UPDATE: Eigen now checks the C++ standard used by CMake. To avoid errors, change the minimum CXX_STANDARD in `lammps/cmake/CMakeLists.txt` from 11 to 17.

The files in `matlab` are for running the sampling in parallel. (Currently, the C++ code can only sample in serial.) Using this code is recommended. To use this code, include the `matlab` directory in your simulation directory. Modify the simulation parameters in `sampling.m` and change the path to the ROBs from your reference potentials and to your samples ROB directory. Then run `sampling.m`.

## Background

The goal of proper orthogonal decomposition is to take a complex system of seemingly random vectors and extract some kind of order from the chaos. This is done by modeling the trajectory of each particle as a function of a spatially-dependent function and a time-dependent coefficient.

$$
\pmb{q}(t) = \pmb{q}(0) + \sum^{\infty}_{k=1} a_k(t) \pmb{\Phi}_k 
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

### fix nvt/rom command
```
fix ID group-ID style_name keyword value ...
```
* ID, group-ID are documented in [fix](https://docs.lammps.org/fix.html) command
* style_name = *nvt/rom* (*npt/rom* and *nph/rom* are not stable)
* one keyword must be appended
```
keyword = model
values = order robfile
  order = model order, which corresponds to the number of columns to be read from the
      reduced-order basis
  robfile = filepath to the input file containing the reduced-order basis
```

The `rom` fixes run a reduced-order simulation in the chosen ensemble. The code reads the reduced-order basis, converts between physical and reduced quantities and calculate time integration in the reduced-order space. The reduced model *should* be compatible with running on multiple processors using MPI.

NOTE: To use the `rom` command with the NVT ensemble, make sure to replace the `fix_nh_rom.cpp` file in the LAMMPS `src` directory. The file includes a few lines of code for parsing the keyword 'model'.

### fix rob command

```
fix ID group-ID rob N order robfile
```
* ID, group-ID are documented in [fix](https://docs.lammps.org/fix.html) command
* rob = style name of this fix command
* N = takes snapshots every N timesteps
* order = model order, which corresponds to the number of columns to be written from the reduced-order basis
* robfile = filepath to the output file for the reduced-order basis
```
keyword = global
  value = integer n that is the number of models (prompts the ROB for the global basis to be calculated after n-1 post_run simulations)
```

The `rob` fix generates a reduced order basis using the POD method. The code assembles an array of snapshots during the run, computes the reduced-order basis and writes the results to a user-designated file. The output file format is space-delimited, with the singular values listed in decreasing order at the bottom of the file.

This command can be used with the LAMMPS [rerun](https://docs.lammps.org/rerun.html) command, as long as the timestep is set to a multiple of the timesteps in the dump files. **To decrease the memory usage during the initial run of the full-simulation, generating the reduced-order basis from a rerun is recommended.** Generating from a rerun also allows the user to create an basis for a simulation run on multiple processors.

### fix rob/stiefel command

```
fix ID group-ID rob/stiefel Nsamples order file1 file2 ... fileM sampleformat keyword value
```
* ID, group-ID are documented in [fix](https://docs.lammps.org/fix.html) command
* Nsamples = number of samples to generate
* order = model order, which corresponds to the number of columns to be written from the reduced-order basis
* file1, file2, ... fileM = ROB files to read, with the last file being the global basis
* sampleformat = format string for generate ROB sample files
* one keyword may be appended
```
keyword = seed
  value = integer seed for the random number generator
keyword = select
  value = none (including this keyword prompts the sampling to also include a "selection.txt" file with the sample-based potentials)
```

The `rob/stiefel` command generates ROB samples based on a few input ROB files. The code projects the bases onto the tangent space of the Stiefel manifold, and randomly samples from the tangent space. The samples are then retracted back to the Stiefel manifold and printed to ROB files. The method is outlined in [Zhang and Guilleminot 2023](https://doi.org/10.1016/j.cma.2022.115702).

## Author

[UQ-ROM](https://github.com/skounouho/UQ-ROM) was created at Sandia National Labs by Senou Kounouho based on work by Hao Zhang and Johann Guilleminot, PhD, at Duke University.

This work is based on previous libraries:
* [eiquadprog.hpp](http://www.cs.cmu.edu/~bstephe1/eiquadprog.hpp), Copyright (2011) Benjamin Stephens, GPL v2.
* [LAMMPS](http://www.lammps.org/), Copyright (2003) Sandia Corporation, GPL v2.

## References

Zhang H., Guilleminot J., A Riemannian stochastic representation for quantifying model uncertainties in molecular dynamics simulations, Comput. Methods Appl. Mech. Engrg., 403 (2023), Article 115702, [10.1016/j.cma.2022.115702](https://doi.org/10.1016/j.cma.2022.115702)

Thompson A.P., Aktulga H.M., Berger R., Bolintineanu D.S., Brown W.M., Crozier P.S., in’t Veld P.J., Kohlmeyer A., Moore S.G., Nguyen T.D., Shan R., Stevens M.J., Tranchida J., Trott C., Plimpton S.J. LAMMPS - a flexible simulation tool for particle-based materials modeling at the atomic, meso, and continuum scales
Comput. Phys. Comm., 271 (2022), Article 108171, [10.1016/j.cpc.2021.108171](https://doi.org/10.1016/j.cpc.2021.108171)
