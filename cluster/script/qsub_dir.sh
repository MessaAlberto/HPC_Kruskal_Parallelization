#!/bin/bash

#PBS -l select=8:ncpus=64:mem=128gb
#PBS -q short_cpuQ
#PBS -l walltime=06:00:00

module load mpich-3.2

DIR=${PBS_O_WORKDIR}

mpiexec -n 512 "${DIR}/bin/kruskal" $V