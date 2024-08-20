#!/bin/bash

FILE=a.out
if test -f "$FILE"; then
    rm a.out
fi

gfortran -O3 -ffast-math -funroll-loops params.f90 functions.f90 propagator.f90 main.f90 -fopenmp
./a.out
rm a.out
rm *.mod