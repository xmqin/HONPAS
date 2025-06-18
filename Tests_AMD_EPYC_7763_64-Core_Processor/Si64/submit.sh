#!/bin/bash
#SBATCH -J dghf
#SBATCH -p normal
#SBATCH -n 128
#SBATCH -N 2
#SBATCH --time=240:00:00
#SBATCH -o job-%J.out
#SBATCH -e job-%J.err
#SBATCH --qos=cpu
#SBATCH --exclusive

# load the environment
#module purge
#module load gcc/9.3.0
#source /public/software/compiler/intel/oneapi2022/setvars.sh intel64
#export I_MPI_DEBUG=5
#export I_MPI_PIN_PROCESSOR_LIST=all
export I_MPI_FABRICS=shm:ofi
export FI_PROVIDER=verbs
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/public/home/xmqin/MathLibs/libxc/6.2.2/intel2020/lib:/public/home/xmqin/MathLibs/fftw/3.3.10/intel2020/lib
#export I_MPI_PMI_LIBRARY=/opt/gridview/slurm/lib/libpmi2.so
#export PATH=/public/software/apps/vasp/5.4.4/hpcx-2.4.1-intel2017:${PATH}
#APP=/public/home/xmqin/git/pwdft_install/examples/pwdft
APP=/public/home/xmqin/HONPAS/honpas_v2.0/Obj/honpas
# run vasp on 2 nodes with 16 cores
mpirun $APP  <Sibulk.fdf > Sibulk.out

