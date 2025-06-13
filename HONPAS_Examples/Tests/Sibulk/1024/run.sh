#!/bin/bash
#SBATCH -p G1Part_sce
#SBATCH -N 5
#SBATCH --ntasks-per-node=48
#SBATCH -n 240
#SBATCH --exclusive
#SBATCH -J honpas
#SBATCH -t 240:00:00

unset I_MPI_PMI_LIBRARY
export I_MPI_JOB_RESPECT_PROCESS_PLACEMENT=0

ulimit -s unlimited
ulimit -l unlimited
module purdge
module add intel/18
module add mpi/intel/18

HONPAS=/es01/paratera/sce4927/xmqin/intel2018/honpas_2.0-beta/Obj/honpas
mpirun -n 240 $HONPAS <Sibulk.fdf > Sibulk.out


