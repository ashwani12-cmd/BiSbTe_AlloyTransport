#!/bin/bash
#SBATCH --job-name=gpumd
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=gpumd.out
#SBATCH --error=gpumd.err

############################
# Environment
############################
module purge
module load gcc
module load cuda/12.8.0-gcc-14.2.0
#module load cuda/12.0

export CUDA_DEVICE_ORDER=PCI_BUS_ID
# DO NOT hardcode CUDA_VISIBLE_DEVICES unless debugging
# export CUDA_VISIBLE_DEVICES=0

echo "=============================================="
echo " Job start time : $(date)"
echo " Hostname       : $(hostname)"
echo " SLURM Job ID   : $SLURM_JOB_ID"
echo "=============================================="

############################
# CUDA / software info
############################
echo "=== CUDA SOFTWARE STACK ==="
which nvcc
nvcc --version
nvidia-smi
echo

############################
# GPU TOPOLOGY
############################
echo "=== GPU TOPOLOGY ==="
nvidia-smi topo -m
echo

############################
# MIG STATUS (CRITICAL)
############################
echo "=== MIG STATUS ==="
nvidia-smi -L
nvidia-smi -q | sed -n '/MIG Mode/,+10p'
echo

############################
# DETAILED GPU INFO
############################
echo "=== DETAILED GPU INFO ==="
nvidia-smi -q -d \
MEMORY,UTILIZATION,POWER,CLOCK,COMPUTE,ECC,PCI,PERFORMANCE
echo

############################
# CLOCKS & POWER LIMITS
############################
echo "=== CLOCKS & POWER ==="
nvidia-smi --query-gpu=index,name,pstate,clocks.sm,clocks.mem,clocks.gr,clocks.max.sm,clocks.max.mem,power.limit,power.draw \
--format=csv
echo

############################
# ENVIRONMENT
############################
echo "=== ENVIRONMENT VARIABLES ==="
echo "CUDA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES"
echo "LD_LIBRARY_PATH=$LD_LIBRARY_PATH"
echo "PATH=$PATH"
echo

############################
# START GPU MONITOR (background)
############################
echo "=== Starting GPU monitor ==="
nvidia-smi dmon -s pucm -d 5 > gpu_monitor.log &
DMON_PID=$!
sleep 2

############################
# RUN GPUMD
############################
cd $SLURM_SUBMIT_DIR

# NEW (Your v4.4 binary):
GPUMD_BIN=/home/IITB/multiscale-mechanics/ashwani12/GPUMD/src/gpumd

echo "=== Starting GPUMD ==="
echo "Start time: $(date)"

srun $GPUMD_BIN > out.dat

GPUMD_EXIT_CODE=$?

echo "GPUMD exit code: $GPUMD_EXIT_CODE"
echo "End time: $(date)"

############################
# STOP GPU MONITOR
############################
kill $DMON_PID
sleep 1

############################
# FINAL GPU STATE
############################
echo "=== FINAL GPU STATE ==="
nvidia-smi
echo

echo "=============================================="
echo " Job finished at: $(date)"
echo "=============================================="


