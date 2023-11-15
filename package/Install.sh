# define LAMMPS directory
lammps_dir="/path/to/lammps"
[ -d $lammps_dir ] || exit 1

# replace fix_nh in lammps/src
if [ -e "$lammps_dir/src/fix_nh.cpp" ]; then
	echo "Replacing $lammps_dir/src/fix_nh.cpp..."
	rm "$lammps_dir/src/fix_nh.cpp"
	cp "./fix_nh.cpp" "$lammps_dir/src/fix_nh.cpp"
else
	echo "Error: $lammps_dir/src/fix_nh.cpp not found."
fi

# add UQ-ROM package to lammps/src
if [ -d "./UQ-ROM" ]; then
	echo "Copying UQ-ROM into $lammps_dir/src..."
	cp -r "./UQ-ROM" "$lammps_dir/src/UQ-ROM"
else
	echo "Error: UQ-ROM not found."
fi

# add UQ-ROM.cmake to lammps/cmake
if [ -e "./cmake/UQ-ROM.cmake" ]; then
	echo "Copying UQ-ROM into $lammps_dir/cmake/Modules/Packages..."
	cp "./cmake/UQ-ROM.cmake" "$lammps_dir/cmake/Modules/Packages/UQ-ROM.cmake"
else
	echo "Error: UQ-ROM.cmake not found"
fi

# add uqrom to lib lammps/lib
if [ -d "lib/uqrom" ]; then
	echo "Copying uqrom into $lammps_dir/lib..."
	cp -r "./lib/uqrom" "$lammps_dir/lib/uqrom"
else
	echo "Error: lib/uqrom not found"
fi

# modify CMakeLists
if [ $(grep -c "UQ-ROM" "$lammps_dir/cmake/CMakeLists.txt") -le 1 ]; then
	echo "Modifying CMakeLists.txt..."
	sed -i 's/UEF/UEF\
  UQ-ROM/' "$lammps_dir/cmake/CMakeLists.txt"
	sed -i "s/MACHDYN VTK/MACHDYN UQ-ROM VTK/" "$lammps_dir/cmake/CMakeLists.txt"
else
	echo "$lammps_dir/cmake/CMakeLists.txt already modified."
fi

# modify Makefile
if [ $(grep -c "uq-rom" "$lammps_dir/src/Makefile") -eq 0 ]; then
	echo "Modifying Makefile..."
	sed -i 's/vtk /vtk \\\
    uq-rom /' "$lammps_dir/src/Makefile"
else
	echo "$lammps_dir/src/Makefile already modified."
fi