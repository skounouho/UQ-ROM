#!/usr/bin/env bash

# define LAMMPS directory
lammps_dir="/home/senou/research/mylammps"

# copy fix files to lammps/src
for f in src/*
do
	echo "Removing $f from lammps/src..."
	rm "$lammps_dir/$f"
done
echo "Fix files deleted"

# restoring fix_nh.cpp
cp "src/fix_nh.cpp" "$lammps_dir/src/fix_nh.cpp"
echo "Restored fix_nh.cpp"