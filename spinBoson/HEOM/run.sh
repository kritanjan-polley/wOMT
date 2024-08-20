#!/bin/bash

FILE=a.out
if test -f "$FILE"; then
    rm a.out
fi

# gfortran -fcheck=all -fbacktrace -Wall -Wextra heom4.f90 spectra_heom4.f90 -fopenmp
gfortran -O3 -ffast-math -funroll-loops heom4.f90 spectra_heom4.f90 -fopenmp
./a.out 
rm a.out *.mod
