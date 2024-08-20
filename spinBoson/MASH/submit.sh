#!/bin/bash
FILE=a.out
if test -f "$FILE"; then
    rm a.out
fi

gfortran -O3 -ffast-math -funroll-loops params.f90 functions.f90 force.f90 modMultiMash.f90 main.f90 -framework Accelerate
./a.out $RANDOM
rm a.out
rm *.mod