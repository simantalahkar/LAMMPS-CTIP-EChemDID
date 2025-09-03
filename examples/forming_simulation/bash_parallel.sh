#!/bin/bash
#SBATCH -J memrstr
#SBATCH --nodes=6
#SBATCH --ntasks-per-node=128
#SBATCH --partition=rome
#SBATCH --mail-user=s.lahkar@tue.nl
#SBATCH --mail-type=BEGIN,END
#SBATCH --time=72:00:00
 
# Load modules for MPI and other parallel libraries
module load 2022
module load libtirpc/1.3.2-GCCcore-11.3.0
module load MPICH/4.0.2-GCC-11.3.0

srun -N 6 $HOME/opt/lammps2018_modified/bin/lmp_mpi -in in.lammps
