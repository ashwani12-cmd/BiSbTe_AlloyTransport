#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=48
#SBATCH --exclusive
#SBATCH --time=01:00:00
#SBATCH --partition=small
#SBATCH --output=C11-%A_%03a.out
#SBATCH --error=C11-%A_%03a.err
#SBATCH --job-name=C11_+0.002
#SBATCH --mem=800

# Load environment
module load openmpi/4.1.4
module load nvhpc/24.11
module load intel-oneapi-mkl/2024.2.2-oneapi-2025.0.1-5u4sz3m

# Set OpenMP threads
export OMP_NUM_THREADS=1

# Path to correct mpirun (HPC-X inside nvhpc)
MPIRUN_PATH="/home/apps/hpc_sdk/Linux_x86_64/24.11/comm_libs/12.6/hpcx/hpcx-2.20/ompi/bin/mpirun"

# Path to your pw.x binary
PWX_BIN="/home/IITB/multiscale-mechanics/amit.k.singh/software/builds/qe-7.1/build/bin/pw.x"

# Move to job submission directory
cd $SLURM_SUBMIT_DIR

# Run Quantum ESPRESSO
time ${MPIRUN_PATH} -np "$SLURM_NTASKS" ${PWX_BIN} < espresso.pwi > espresso.pwo
