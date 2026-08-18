#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=48
#SBATCH --exclusive
#SBATCH --time=01:00:00
#SBATCH --mem=24G
#SBATCH --job-name=lammps
#SBATCH --partition=debug
#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out

# Runtime linking and compile-time flags for libevent & hwloc
export LDFLAGS="-L/home/apps/spack/opt/spack/linux-almalinux8-cascadelake/gcc-14.2.0/libevent-2.1.12-2fqz3dhvnca3g6neeavbywpxp5z7fy4v/lib \
                -L/home/apps/spack/opt/spack/linux-almalinux8-cascadelake/gcc-14.2.0/hwloc-2.11.1-ma4yhksnj74h67tta3vtf4eqloilx6ne/lib"
export CPPFLAGS="-I/home/apps/spack/opt/spack/linux-almalinux8-cascadelake/gcc-14.2.0/libevent-2.1.12-2fqz3dhvnca3g6neeavbywpxp5z7fy4v/include \
                 -I/home/apps/spack/opt/spack/linux-almalinux8-cascadelake/gcc-14.2.0/hwloc-2.11.1-ma4yhksnj74h67tta3vtf4eqloilx6ne/include"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/home/apps/spack/opt/spack/linux-almalinux8-cascadelake/gcc-14.2.0/libevent-2.1.12-2fqz3dhvnca3g6neeavbywpxp5z7fy4v/lib:/home/apps/spack/opt/spack/linux-almalinux8-cascadelake/gcc-14.2.0/hwloc-2.11.1-ma4yhksnj74h67tta3vtf4eqloilx6ne/lib"

# Load compiler and runtime modules used during build
module purge
module load glibc/2.28-gcc-8.5.0-bfcbr6w
module load gcc-runtime/8.5.0-gcc-8.5.0-cqzkmjo
module load zlib-ng/2.2.1-gcc-8.5.0-orhazg6
module load binutils/2.43.1-gcc-8.5.0-kdsptcv
module load zstd/1.5.6-gcc-8.5.0-fg5aduu
module load gmp/6.3.0-gcc-8.5.0-4gikgah
module load mpfr/4.2.1-gcc-8.5.0-wvjtu4c
module load mpc/1.3.1-gcc-8.5.0-qiu6dv2
module load gcc/14.2.0-gcc-8.5.0-777kyuf
module load openmpi/4.1.4

# Backup loaded modules for reproducibility
module list > modules_used_${SLURM_JOB_ID}.txt

# Ensure single-threaded OpenMP per MPI rank
export OMP_NUM_THREADS=1

# Avoid stack overflows (MLIP and deep models can cause this)
ulimit -s unlimited

# Move to submission directory
cd $SLURM_SUBMIT_DIR

# Informative log
echo "Running on $SLURM_NTASKS tasks across $SLURM_NNODES nodes"
echo "Starting LAMMPS job at $(date)"

# Run LAMMPS
time mpirun --mca btl ^openib -np $SLURM_NTASKS /home/IITB/multiscale-mechanics/ashwani12/software/builds/NEP_CPU/lammps-24Mar2022/src/lmp_mpi -in in.elastic > out.dat


echo "Job finished at $(date)"

