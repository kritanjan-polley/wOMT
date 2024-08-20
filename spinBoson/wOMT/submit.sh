#!/bin/bash

FILE=a.out
if test -f "$FILE"; then
    rm a.out
fi

gfortran -O3 -ffast-math -funroll-loops params.f95 functions.f95 propagator.f95 main.f95 -fopenmp
./a.out
rm a.out
rm *.mod
