#!/bin/bash
#PBS -N gpumd_job
#PBS -q gpu
#PBS -l select=1:ncpus=8:ngpus=1:mem=64gb
#PBS -l walltime=96:00:00
#PBS -o output.log
#PBS -e error.log

echo "================= JOB INFO ================="
echo "Host: $(hostname)"
echo "Job ID: $PBS_JOBID"
echo "Nodefile:"
cat $PBS_NODEFILE
echo "Working Directory: $PBS_O_WORKDIR"
echo "============================================"

cd "$PBS_O_WORKDIR"

# ---------------- Environment Setup ----------------
if [ -f /apps/spack/v0.23/share/spack/setup-env.sh ]; then
    . /apps/spack/v0.23/share/spack/setup-env.sh
    spack load cuda@12.6.2%gcc@11.4.1 arch=linux-rhel9-icelake
fi

# -------- USE ONLY 3 GPUs (VERY IMPORTANT) ---------
export CUDA_VISIBLE_DEVICES=0
echo "CUDA_VISIBLE_DEVICES = $CUDA_VISIBLE_DEVICES"

echo "------------- GPU STATUS BEFORE RUN -------------"
nvidia-smi
echo "--------------------------------------------------"

# ----------------- Run GPUMD ------------------
echo "Starting GPUMD run..."
/home/multiscale-mechanics/ashwani12/software/builds/GPUMD/src/gpumd > out.dat

echo "================ JOB FINISHED ================="

