#!/bin/bash
#SBATCH -p G1Part_sce
#SBATCH -N 8
#SBATCH --ntasks-per-node=48
#SBATCH -n 384
#SBATCH --exclusive
#SBATCH -J honpas
#SBATCH -t 240:00:00

unset I_MPI_PMI_LIBRARY
export I_MPI_JOB_RESPECT_PROCESS_PLACEMENT=0

ulimit -s unlimited
ulimit -l unlimited
module purdge
source /es01/paratera/parasoft/oneAPI/2022.1/setvars.sh
HONPAS=/es01/paratera/sce4927/xmqin/intel2022/honpas_2.0-beta/Obj/honpas
mpirun -n 384 $HONPAS <tio2.fdf > tio2.out


