# Examples

Example inputs for simulations with the integrated CTIP–EChemDID package, and supporting scripts.

- `forming_simulation/`: a LAMMPS input script (`in.lammps`, written for LAMMPS 22 Aug 2018) for a forming simulation of a TaOx/HfO2 device under applied bias, with its relaxed starting structure, the CTIP force-field files, a SLURM job script and a Python script that tracks filament evolution in the output trajectories with [LAMMPSKit](https://github.com/simantalahkar/lammpskit).
- `forcefield_coefficient_generation/`: the script used to generate the force-field coefficient tables (legacy research code) and the generated tables.
