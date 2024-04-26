#!/usr/bin/env bash

# define LAMMPS directory
lammps_dir="/home/senou/research/mylammps"
eigen_dir="/home/senou/research/eigen"

# copy fix files to lammps/src
for f in src/*
do
	echo "Linking $f ..."
	ln $f "$lammps_dir/$f"
done
echo "Fix files copied"

# symlink eigen
echo "Linking Eigen ..."
ln -s "$eigen_dir" "$lammps_dir/src/eigen"
echo "Eigen library linked"