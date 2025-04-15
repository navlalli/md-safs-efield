[![CC BY 4.0][cc-by-shield]][cc-by]

# Modelling combustion in external electric fields with LAMMPS

This repository contains the [LAMMPS](https://www.lammps.org/) input scripts and force fields required to reproduce the results presented in *On the modelling of hydrocarbon combustion in external electric fields with reactive molecular dynamics*. 

LAMMPS from 2 April 2025 or later (see the [build](build) directory for package information) is required to run the input scripts.

The following table describes the contents of this repository:

| Directory | Description |
| :--- | :--- |
| [build](./build) | Build information |
| [reaxff](./reaxff) | Force fields and Gaussian exponents |
| [q_carbondioxide_water](./q_carbondioxide_water) | Charge analysis of one water molecule and one carbon dioxide molecule |
| [q_water_two](./q_water_two) | Charge analysis of two water molecules |
| [q_pmh_one](./q_pmh_one) | Charge analysis of one pentamethylheptane molecule |
| [q_pmh_oxygen](./q_pmh_oxygen) | Charge analysis of pentamethylheptane in oxygen |
| [pmh_combustion](./pmh_combustion) | Equilibration and combustion of pentamethylheptane in oxygen |
| [post_process](./post_process) | Python post-processing script to analyze kinetic energies |


## License

This work is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by].

[![CC BY 4.0][cc-by-image]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg
