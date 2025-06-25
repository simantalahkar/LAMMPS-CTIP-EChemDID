# CTIP–EChemDID Integrated Source Files

This directory contains the modified source code files for the integrated CTIP–EChemDID package for LAMMPS.

## Overview

- **CTIP (Charge-Transfer Ionic Potential)**: Implements dynamic charge transfer using QEq with the Streitz-Mintmire pair style. Originally based on work by Zhou and Wadley.
- **EChemDID (Electrochemical Dynamics with Implicit Degrees of Freedom)**: Adds fictitious dynamics for electrochemical potential equilibration within metallic regions, dynamically affecting atomic electronegativities.

This integration allows both models to be used concurrently in LAMMPS simulations. The code has been adapted and tested to work for the 22Aug18 stable version of LAMMPS, and dependencies between packages have been resolved.

## Modifications by Simanta Lahkar (2025)

- Updated CTIP for compatibility with the 2018 LAMMPS version.
- Modified CTIP and EChemDID components to interoperate with each other: CTIP qeq uses the local potential from EChemDID to equilibrate the charges.

## Installation

1. Enable required packages (depending on your LAMMPS setup):

   ```bash
   make yes-qeq
   make yes-user-misc
   ```

   Or with CMake:

   ```bash
   cmake ../cmake -D PKG_QEQ=on -D PKG_USER_MISC=on
   make -j4
   ```

2. Copy the files from the CTIP and EChemDID folders into the `src/` directory of LAMMPS:

   ```bash
   cp -r CTIP-EChemDID-Integrated-22Aug18/CTIP/* /$user-directory-path/lammps/src/
   cp -r CTIP-EChemDID-Integrated-22Aug18/EChemDID/* /$user-directory-path/lammps/src/
   ```

3. Build LAMMPS as usual.

## Citation Prompt

During initialization of a LAMMPS run using this package, a message will be printed asking the user to cite the relevant papers listed below.

## Citation Instructions

Please cite the following works if you use this code:

- **This Integration**:
  - S. Lahkar et al., Decoupling Local Electrostatic Potential and Temperature-Driven Atomistic Forming Mechanisms in TaOx/HfO2-Based ReRAMs using Reactive Molecular Dynamics Simulations, arXiv:2505.24468 or https://doi.org/10.48550/arXiv.2505.24468

- **CTIP**:
  - Chem. Mater. 2017, 29, 8, 3603–3614
  - Chem. Mater. 2019, 31, 9, 3089–3102

- **EChemDID**:
  - Onofrio & Strachan, *J. Chem. Phys.* 143, 054109 (2015) — https://doi.org/10.1063/1.4927562  
  - Onofrio et al., *Nature Materials* 14, 4 (2015) — https://doi.org/10.1038/NMAT4221

## Contact

For issues, suggestions, or collaboration requests, please contact:

**Simanta Lahkar**  
Email: s.lahkar@tue.nl, simantalahkar@hotmail.com 
GitHub: https://github.com/simantalahkar

## License

This work is distributed under the GNU General Public License (GPL). See the `LICENSE` file in the root directory for full terms.
