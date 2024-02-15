# define LAMMPS directory
lammps_dir="/path/to/lammps"
[ -d $lammps_dir ] || exit 1

# remove UQ-ROM package from lammps/src
if [ -d "./UQ-ROM" ]; then
	echo "Removing UQ-ROM from $lammps_dir/src..."
	rm -r "$lammps_dir/src/UQ-ROM"
else
	echo "Error: UQ-ROM not found."
fi

# remove UQ-ROM.cmake from lammps/cmake
if [ -e "./cmake/UQ-ROM.cmake" ]; then
	echo "Removing UQ-ROM from $lammps_dir/cmake/Modules/Packages..."
	rm "$lammps_dir/cmake/Modules/Packages/UQ-ROM.cmake"
else
	echo "Error: UQ-ROM.cmake not found"
fi

# remove uqrom from lib lammps/lib
if [ -d "lib/uqrom" ]; then
	echo "Removing uqrom from $lammps_dir/lib..."
	rm -r "$lammps_dir/lib/uqrom"
else
	echo "Error: lib/uqrom not found"
fi

# restore CMakeLists
if [ $(grep -c "UQ-ROM" "$lammps_dir/cmake/CMakeLists.txt") -ge 1 ]; then
	echo "Restoring CMakeLists.txt..."
	sed -i 's/UEF\
  UQ-ROM/UEF/' "$lammps_dir/cmake/CMakeLists.txt"
	sed -i "s/MACHDYN UQ-ROM VTK/MACHDYN VTK/" "$lammps_dir/cmake/CMakeLists.txt"
else
	echo "$lammps_dir/cmake/CMakeLists.txt already restored."
fi

# restore Makefile
if [ $(grep -c "uq-rom" "$lammps_dir/src/Makefile") -ge 1 ]; then
	echo "Restoring Makefile..."
	sed -i 's/vtk \\\
    uq-rom /vtk /' "$lammps_dir/src/Makefile"
else
	echo "$lammps_dir/src/Makefile already restored."
fi
