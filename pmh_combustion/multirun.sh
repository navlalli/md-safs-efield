#!/bin/bash

# Set method to compute charges
qcalc="qeq" # qeq, qeqr or qtpie (qeq only used when no external electric field)

# Set method for thermostatting or nve
tmode="nhoover" # nhoover or csvr or nve

# Set efield
ef="0.0" # z-cpt of electric field (V/Angstrom)

# Create case directory
case_index=1
case_name="${qcalc}_run_${tmode}_temp2000/ef_${ef}/run_${case_index}"
mkdir -p ${case_name}
cd ${case_name}
echo "Running ${case_name}"

# Run LAMMPS
mpirun -n 12 ~/lammps-2Apr2025/build/lmp -nocite -log log.run \
	-var qcalc ${qcalc} \
	-var tmode ${tmode} \
	-var case_ind ${case_index} \
	-var ef ${ef} \
	-in ../../../in.multirun
