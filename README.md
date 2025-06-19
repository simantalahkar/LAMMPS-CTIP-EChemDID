# LAMMPS-CTIP-EChemDID
This repository contains an integrated and updated implementation of two LAMMPS packages:

- **CTIP** (Charge-Transfer Ionic Potential) - CTIP modifications to the Streitz-Mintmire potential and QEq-based charge equilibration.
- **EChemDID** (Electrochemical Dynamics with Implicit Degrees of Freedom) - description of local electrochemical potentials and their equilibration through metallic atoms with fictitious dynamics.

The purpose of this integration is to enable combined simulations involving both electrochemical potential equilibration and dynamic charge-transfer interactions in metal oxide materials. The packages have been adapted for compatibility with a newer version of LAMMPS and modified to work together.

## What's Included

- Combined and modified source files for both CTIP and EChemDID in the `CTIP-ECHEMDID-Integrated-22Aug18/` folder.
- Updated code for compatibility with modern LAMMPS builds.
- Instructions for compilation and usage.

## Author of Integration

**Simanta Lahkar**, 2025  
Eindhoven University of Technology
Contact: s.lahkar@tue.nl, simantalahkar@hotmail.com 

## Citing This Work

If you use this integrated package for simulation or publication, please cite the following:

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

This project is licensed under the GNU General Public License (GPL), consistent with the original CTIP and EChemDID repositories. See the `LICENSE` file for full details.

