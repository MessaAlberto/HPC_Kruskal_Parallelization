#!/bin/bash

#PBS -l select=2:ncpus=64:mem=128gb
#PBS -q short_cpuQ
#PBS -l walltime=00:10:00

mpiexec -n 128 ./cluster/kruskal $V