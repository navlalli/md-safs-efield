#!/bin/bash

qcalc="qeqr"  # qeqa, qeqr or qtpie
scale="1"  # Investigated scale=1 and scale=8
mkdir -p "${qcalc}_scale${scale}"
cd "${qcalc}_scale${scale}"

mpirun -np 1 ~/lammps-2Apr2025/build/lmp -var qcalc ${qcalc} -var scale ${scale} \
       -in ../in.pmh -log log.run
