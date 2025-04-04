#!/bin/bash

qcalc="qtpie"  # qeqa, qeqr or qtpie
mkdir -p ${qcalc}
cd ${qcalc}

mpirun -np 12 ~/lammps-2Apr2025/build/lmp -var qcalc ${qcalc} -in ../in.pmh_oxygen -log log.run
