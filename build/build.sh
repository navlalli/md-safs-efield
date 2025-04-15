#!/bin/bash

# Move this file to main lammps directory, such as ~/lammps-2Apr2025/ , and run to build lammps using cmake
mkdir build
cd build
cmake ../cmake
cmake -C ../cmake/presets/basic.cmake .
cmake -D PKG_REAXFF=ON .
cmake -D PKG_EXTRA-FIX=ON .

cmake --build .
