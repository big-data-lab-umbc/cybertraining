#!/bin/tcsh

#f2py -m add add.f
#f2py3 -c -m calc_dist calc_dist.f90
f2py3 -c --fcompiler=gfortran --f90flags='-fopenmp' -lgomp -m k_means_mod k-means_mod.f90
