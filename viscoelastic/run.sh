#!/bin/bash
 We="1.0"
 Ec="0.1"
 De="1.0"

qcc -fopenmp -Wall -O2 -disable-dimensions bounceVE.c -o bounceVE -lm
export OMP_NUM_THREADS=8
./bounceVE $We $Ec $De