#!/bin/bash

qcalc="qeqa"  # qeqa, qeqr or qtpie
mkdir -p ${qcalc}
cd ${qcalc}

mpirun -np 1 ~/lammps-2Apr2025/build/lmp -var qcalc ${qcalc} -in ../in.water -log log.run
