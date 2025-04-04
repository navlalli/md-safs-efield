#!/bin/bash

qcalc="qtpie"  # qeqa, qeqr or qtpie
mkdir -p ${qcalc}
cd ${qcalc}

dzs=(0.0 25.0 49.5)
ezs=(0.05)
for dz in ${dzs[*]}; do
    for ez in ${ezs[*]}; do
	mpirun -np 1 ~/lammps-2Apr2025/build/lmp -in ../in.co2_h2o -var dz ${dz} \
	       -var ez ${ez} -var qcalc ${qcalc} -log log.dz${dz}_ez${ez}.${mode}
    done
done
