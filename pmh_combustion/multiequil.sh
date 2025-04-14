#!/bin/bash

# Set method for computing charges
qcalc="qeq" # qeq or qtpie

# Set method for thermostatting
tmode="nhoover" # nhoover or csvr

# Create case directory
temp="1500" # 1500 or 2000 K
case_index=1
case_name="${qcalc}_equil_${tmode}_temp${temp}/equil_${case_index}"
mkdir -p ${case_name}
cd ${case_name}
echo "Running ${case_name}"

# Run LAMMPS
mpirun -n 12 ~/lammps-2Apr2025/build/lmp -nocite -log log.equil \
	-var qcalc ${qcalc} \
	-var tmode ${tmode} \
	-var case_ind ${case_index} \
	-in ../../in.multiequil_temp${temp}
