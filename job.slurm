#!/bin/bash

#SBATCH -N 2
#SBATCH -n 2
#SBATCH -c 16
#SBATCH -t 00:20:00

module load openmpi/1.8.5
module load intel/2015
mpic++ main.cc -std=c++11 -Wall -openmp -O3 -xHost

# try running using different commands

export OMP_NUM_THREADS=1
echo " "
echo "Running with OMP_NUM_THREADS=1"
echo "------------------------------"
echo " "
echo "Running normally..."
time ./a.out
echo " "
echo "Running using srun..."
time srun ./a.out
echo " "
echo "Running using mpirun..."
time mpirun ./a.out
echo " "
echo "Running using mpiexec..."
time mpiexec ./a.out

export OMP_NUM_THREADS=2
echo " "
echo "Running with OMP_NUM_THREADS=2"
echo "------------------------------"
echo " "
echo "Running normally..."
time ./a.out
echo " "
echo "Running using srun..."
time srun ./a.out
echo " "
echo "Running using mpirun..."
time mpirun --map-by ppr:2:socket ./a.out
echo " "
echo "Running using mpiexec..."
time mpiexec ./a.out

export OMP_NUM_THREADS=4
echo " "
echo "Running with OMP_NUM_THREADS=4"
echo "------------------------------"
echo " "
echo "Running normally..."
time ./a.out
echo " "
echo "Running using srun..."
time srun ./a.out
echo " "
echo "Running using mpirun..."
time mpirun ./a.out
echo " "
echo "Running using mpiexec..."
time mpiexec ./a.out

export OMP_NUM_THREADS=8
echo " "
echo "Running with OMP_NUM_THREADS=8"
echo "------------------------------"
echo " "
echo "Running normally..."
time ./a.out
echo " "
echo "Running using srun..."
time srun ./a.out
echo " "
echo "Running using mpirun..."
time mpirun ./a.out
echo " "
echo "Running using mpiexec..."
time mpiexec ./a.out
