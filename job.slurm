#!/bin/bash

#SBATCH -N 2
#SBATCH -c 11
#SBATCH -t 00:20:00

module load openmpi-gnu/1.8.5
mpic++ main.cc --std=c++11 -Wall -fopenmp -O3

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
time mpirun ./a.out
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
