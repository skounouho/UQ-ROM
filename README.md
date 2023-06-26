# Reduced-Order Modeling for LAMMPS Using Proper Orthogonal Decomposition

This collection of fixes allow for the creation of reduced order models in LAMMPS using POD linear subspace methods. For a primer on the POD method, see [Weiss 2019](https://doi.org/10.2514/6.2019-3333).

## Installation

Download the files and insert them in the `src` folder in your LAMMPS installation. You will need to install [Eigen](https://eigen.tuxfamily.org) to use `fix rob`. See the [LAMMPS documentation](https://eigen.tuxfamily.org) for how to compile LAMMPS with Eigen.

## Usage

### fix nve/rom command

```
fix ID group-ID nve/rom order robfile
```
* ID, group-ID are documented in [fix](https://docs.lammps.org/fix.html) command
* nve/rom = style name of this fix command
* order = model order, which corresponds to the number of columns to be read from the reduced-order basis
* robfile = filepath to the reduced-order basis

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
  order = model order, which corresponds to the number of columns to be read from the reduced-order basis
  robfile = filepath to the reduced-order basis
```

The `rom` fixes run a reduced-order simulation in the chosen ensemble. The code reads the reduced-order basis, converts between physical and reduced quantities and calculate time integration in the reduced-order space.

### fix rob command
