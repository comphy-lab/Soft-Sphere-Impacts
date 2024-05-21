#!/bin/bash
 We="20"

qcc -fopenmp -Wall -O2 -disable-dimensions bounceVE.c -o bounceVE -lm
export OMP_NUM_THREADS=8
./bounceVE $We