#!/bin/bash

#PBS -l select=2:ncpus=64:mem=128gb
#PBS -q short_cpuQ
#PBS -l walltime=00:10:00

module load mpich-3.2

DIR=${PBS_O_WORKDIR}

mpiexec -n 128 "${PBS_O_WORKDIR}/kruskal" $V